#ifndef NEURONS_HPP
#include <mecacell/mecacell.h>
#include <mecacell/utilities/utils.h>
#include <mecacell/utilities/obj3D.hpp>
#include <cmath>

enum NType { NTinput, NThidden, NToutput };

struct Neuron {
	double nt;
	double input = 0.0;
	double dopamine = 0.0;
	NType type = NThidden;
	int t_fired = 0;
	MecaCell::Vec position;
  double decay = 1.0;
  double delay = 0.0;

	Neuron(double n, NType t, MecaCell::Vec pos) : nt(n), type(t), position(pos) {}

	void step(int t, double vt, double vr) {
		nt += input;
		if (nt >= vt) {
			nt = vr;
			t_fired = t;
		}
		input = 0.0;
	}
};

std::ostream& operator<<(std::ostream& out, const Neuron& n) {
	return out << n.nt << ", " << n.t_fired << ", " << n.type << ", " << n.input << ", " <<
		n.position.coords[0] << ", " << n.position.coords[1] << ", " << n.position.coords[2];
}

struct SNN {
	vector<Neuron> neurons;
	vector<vector<double>> synapses;
	double input_signal;
	double vt;
	double vr;
	double aplus;
	double aminus;

	SNN() {}
	SNN(double i, double t, double r, double p, double m) :
		input_signal(i), vt(t), vr(r), aplus(p), aminus(m) {}

	void fire(int t, double reward_signal) {
		for (size_t i=0; i<neurons.size(); i++) {
			neurons[i].step(t, vt, vr);
		}

		double total_delta = 0.0;
		for (size_t i=0; i<neurons.size(); i++) {
			if (neurons[i].t_fired == t) {
				for (size_t j=0; j<neurons.size(); j++) {
					double s = synapses[j][i];
					if (s > 0.0) neurons[j].input += s;
				}
			}
		}

		for (size_t i=0; i<neurons.size(); i++) {
			if (neurons[i].t_fired == t) {
				for (size_t j=0; j<neurons.size(); j++) {
					double s = synapses[j][i];
					if (s > 0.0) {
						double coeff = (neurons[i].dopamine + neurons[j].dopamine)/2.0;
						double diff = 0.0;
						if (std::abs(neurons[i].t_fired - neurons[j].t_fired) < 10) {
							if (neurons[i].t_fired>neurons[j].t_fired) {
								diff += coeff * aplus*s*(1-s);
							} else if (neurons[i].t_fired<neurons[j].t_fired) {
								diff -= coeff * aminus*s*(1-s);
							}
							synapses[j][i] += diff;
							total_delta += std::abs(diff);
						}
					}
				}
			}
		}
		MecaCell::logger<MecaCell::DBG>("STDP delta :: ", total_delta);

		for (size_t i=0; i<neurons.size(); i++) {
			if (neurons[i].type == NTinput) neurons[i].input = input_signal;
			neurons[i].input += reward_signal;
		}
	}

	void dopamine_release(vector<double> &reward, double da, double dd, double dpa) {
		for (size_t i=0; i<neurons.size(); i++) {
			if (neurons[i].type == NToutput) {
				int x0 = std::floor(neurons[i].delay);
				if (x0 == 10) x0--;
				double interp_reward = reward[x0] + (neurons[i].delay - x0) * (reward[x0+1] - reward[x0]);
				neurons[i].dopamine = (1.0-da)*neurons[i].dopamine + da*neurons[i].decay*interp_reward;
			} else {
				neurons[i].dopamine = (1.0-da)*neurons[i].dopamine + da*reward[0];
			}
		}
	}

	void train() {}
};

std::ostream& operator<<(std::ostream& out, const SNN& s) {
	std::ostringstream sn("");
	sn << "params:\t" << s.input_signal << ", " << s.vt << ", " << s.vr << ", " << s.aplus << ", " << s.aminus;
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
