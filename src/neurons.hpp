#ifndef NEURONS_HPP
#include <mecacell/mecacell.h>
#include <mecacell/utilities/utils.h>
#include <mecacell/utilities/obj3D.hpp>
#include <cmath>

enum NType { NTinput, NThidden, NToutput };

struct Neuron {
	double v;
	double u = 0.0;
  double stdp = 0.0;
	double input = 0.0;
	double dopamine = 0.0;
  bool fired = false;
	NType type = NThidden;
	MecaCell::Vec position = MecaCell::Vec(0,0,0);
  double decay = 1.0;
  double delay = 0.0;
  bool inhibitory = false;

	Neuron(double n, NType t, MecaCell::Vec pos) : v(n), type(t), position(pos) {}

	void step(double vt, double vr) {
    fired = false;
    stdp = stdp * 0.95;
		v = v+0.5*((0.04*v+5)*v+140-u+input);
		v = v+0.5*((0.04*v+5)*v+140-u+input);
    u = u+0.2*(0.2*v-u);
		if (v >= vt) {
			v = vr;
      if (inhibitory) {
        u += 2;
      } else {
        u += 8;
      }
      stdp = 0.1;
      fired = true;
		}
		input = 0.0;
	}
};

std::ostream& operator<<(std::ostream& out, const Neuron& n) {
	return out << n.v << ", " << n.u << ", " << n.type << ", " << n.input << ", " <<
		n.position.coords[0] << ", " << n.position.coords[1] << ", " << n.position.coords[2];
}

struct SNN {
	vector<Neuron> neurons;
	vector<vector<double>> synapses;
	vector<vector<double>> synapses_delta;
	double vt;
	double vr;
	double aplus;
	double aminus;

	SNN() {}
	SNN(double t, double r, double p, double m) : vt(t), vr(r), aplus(p), aminus(m) {}

	void fire(std::vector<double> inputs, double reward_signal) {
		// MecaCell::logger<MecaCell::DBG>("stepping");
		for (size_t i=0; i<neurons.size(); i++) {
			neurons[i].step(vt, vr);
		}

		// MecaCell::logger<MecaCell::DBG>("feed forward input");
		for (size_t i=0; i<neurons.size(); i++) {
			if (neurons[i].fired) {
				for (size_t j=0; j<neurons.size(); j++) {
					neurons[j].input += synapses[i][j];
				}
			}
		}

		// MecaCell::logger<MecaCell::DBG>("next input");
		for (size_t i=0; i<neurons.size(); i++) {
      neurons[i].input += inputs[i];
      neurons[i].input += reward_signal*neurons[i].dopamine;
		}
		// MecaCell::logger<MecaCell::DBG>("done firing");
	}

  void train(double da_factor) {
		// MecaCell::logger<MecaCell::DBG>("STDP");
		for (size_t i=0; i<neurons.size(); i++) {
			if (neurons[i].fired) {
				for (size_t j=0; j<neurons.size(); j++) {
					double s = synapses[j][i];
          // pre-synaptic STDP (LTP)
 					if (synapses[j][i] > 0.0) {
            synapses_delta[j][i] += aplus * neurons[j].stdp;
          }
          // post-synaptic STDP (LTD)
					if (synapses[i][j] > 0.0) {
            synapses_delta[i][j] -= aminus * neurons[i].stdp;
          }
        }
      }
    }

		// MecaCell::logger<MecaCell::DBG>("weight update");
		for (size_t i=0; i<neurons.size(); i++) {
      for (size_t j=0; j<neurons.size(); j++) {
        double s = synapses[i][j];
        if (s > 0.0) {
          double sd = synapses_delta[i][j];
          double da_coeff = (neurons[i].dopamine + neurons[j].dopamine)/2.0;
          double sp = s + (0.01*(1.0-da_factor) + da_factor*da_coeff) * sd;
          sp = max(0.00001, min(4.0, sp));
          synapses[i][j] = sp;
        }
      }
    }

		for (size_t i=0; i<neurons.size(); i++) {
      for (size_t j=0; j<neurons.size(); j++) {
        synapses_delta[i][j] *= 0.9;
      }
    }
  }

  double get_da_mean() {
    double da_mean = 0.0;
    for (size_t i=0; i<neurons.size(); i++) {
      da_mean += neurons[i].dopamine;
    }
    return da_mean / neurons.size();
  }

  double get_da_std(double da_mean) {
    double da_std = 0.0;
    for (size_t i=0; i<neurons.size(); i++) {
      da_std += pow(neurons[i].dopamine - da_mean, 2);
    }
    return sqrt(da_std / neurons.size());
  }

  double get_weight_mean() {
    double weight_mean = 0.0;
    int n_synapse = 0;
    for (size_t i=0; i<neurons.size(); i++) {
      if (!neurons[i].inhibitory) {
        for (size_t j=0; j<neurons.size(); j++) {
          if (synapses[i][j] > 0.0) {
            weight_mean += synapses[i][j];
            n_synapse += 1;
          }
        }
      }
    }
    return weight_mean / n_synapse;
  }

  double get_weight_std(double weight_mean) {
    double weight_std = 0.0;
    int n_synapse = 0;
    for (size_t i=0; i<neurons.size(); i++) {
      if (!neurons[i].inhibitory) {
        for (size_t j=0; j<neurons.size(); j++) {
          if (synapses[i][j] > 0.0) {
            weight_std += pow(synapses[i][j] - weight_mean, 2);
            n_synapse += 1;
          }
        }
      }
    }
    return sqrt(weight_std / n_synapse);
  }

	void dopamine_release(vector<double> &reward, double da, double dd, double dpa) {
		for (size_t i=0; i<neurons.size(); i++) {
      neurons[i].dopamine *= (1.0-da);
			if (neurons[i].type == NToutput) {
				int x0 = std::floor(neurons[i].delay);
				if (x0 == 10) x0--;
				double interp_reward = reward[x0] + (neurons[i].delay - x0) * (reward[x0+1] - reward[x0]);
				neurons[i].dopamine += da*neurons[i].decay*interp_reward;
			} else {
				neurons[i].dopamine += da*reward[0];
			}
		}
	}
};

std::ostream& operator<<(std::ostream& out, const SNN& s) {
	std::ostringstream sn("");
	sn << "params:\t" << s.vt << ", " << s.vr << ", " << s.aplus << ", " << s.aminus;
	sn << "\nneurons:\t";
	for (const auto& n : s.neurons) sn << n << ", ";
	for (size_t i=0; i<s.neurons.size(); ++i) {
		sn << "\nweights:\t" << i << ", ";
		for (size_t j=0; j<s.neurons.size(); ++j) {
			sn << s.synapses[i][j] << ", ";
		}
	}
	return out << sn.str();
}

#endif
