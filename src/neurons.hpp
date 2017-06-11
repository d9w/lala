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
	MecaCell::Vec position;
  double decay = 1.0;
  double delay = 0.0;

	Neuron(double n, NType t, MecaCell::Vec pos) : v(n), type(t), position(pos) {}

	void step(double vt, double vr) {
    fired = false;
    stdp = stdp * 0.95;
		v = v+0.5*((0.04*v+5)*v+140-u+input);
		v = v+0.5*((0.04*v+5)*v+140-u+input);
    u = u+0.2*(0.2*v-u);
		if (v >= vt) {
			v = vr;
      u = u+8;
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

	void fire(std::vector<double> inputs, double reward_signal, double da_factor) {
		// MecaCell::logger<MecaCell::DBG>("stepping");
		for (size_t i=0; i<neurons.size(); i++) {
			neurons[i].step(vt, vr);
		}

		// MecaCell::logger<MecaCell::DBG>("feed forward input");
		for (size_t i=0; i<neurons.size(); i++) {
			if (neurons[i].fired) {
				for (size_t j=0; j<neurons.size(); j++) {
					double s = synapses[j][i];
					neurons[j].input += s;
				}
			}
		}

    int total_fired=0;
		// MecaCell::logger<MecaCell::DBG>("STDP");
		for (size_t i=0; i<neurons.size(); i++) {
			if (neurons[i].fired) {
        total_fired += 1;
				for (size_t j=0; j<neurons.size(); j++) {
					double s = synapses[j][i];
          // pre-synaptic STDP (LTP)
 					if (synapses[j][i] > 0.0) {
            synapses_delta[j][i] += aplus * neurons[j].stdp;
          }
          // post-synaptic STDP (LTD)
					if (synapses[i][j] > 0.0) {
            synapses_delta[i][j] += aminus * neurons[i].stdp;
          }
        }
      }
    }

		// MecaCell::logger<MecaCell::DBG>("weight update");
		double total_delta = 0.0;
		double total_da = 0.0;
		double total_sd = 0.0;
		for (size_t i=0; i<neurons.size(); i++) {
      for (size_t j=0; j<neurons.size(); j++) {
        double s = synapses[i][j];
        if (s > 0.0) {
          double sd = synapses_delta[i][j];
          double da_coeff = (neurons[i].dopamine + neurons[j].dopamine)/2.0;
          double sp = s + (0.01*(1.0-da_factor) + da_factor*da_coeff) * sd;
          sp = max(0.00001, min(2.0, sp));
          total_delta += std::abs(s - sp);
          total_da += da_coeff;
          total_sd += sd;
          synapses[i][j] = sp;
        }
      }
    }

		for (size_t i=0; i<neurons.size(); i++) {
      for (size_t j=0; j<neurons.size(); j++) {
        synapses_delta[i][j] *= 0.9;
      }
    }

    double total_input = 0.0;
		// MecaCell::logger<MecaCell::DBG>("next input");
		for (size_t i=0; i<neurons.size(); i++) {
			if (neurons[i].type == NTinput || neurons[i].type == NThidden) {
        neurons[i].input = inputs[i];
        total_input += inputs[i];
      }
			neurons[i].input += reward_signal*neurons[i].dopamine;
		}

    total_delta/=neurons.size();
    total_da/=neurons.size();
    total_sd/=neurons.size();
    total_input/=neurons.size();
		MecaCell::logger<MecaCell::DBG>("STDP :: ", total_delta, " DA :: ", total_da,
                                    " SD :: ", total_sd, " fired :: ", total_fired,
                                    " input :: ", total_input);


		// MecaCell::logger<MecaCell::DBG>("done firing");
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
