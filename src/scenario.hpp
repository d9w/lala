#ifndef SCENARIO_HPP
#include <mecacell/mecacell.h>
#include <mecacell/utilities/utils.h>
#include <mecacell/plugins/simplifiedfluidplugin.hpp>
#include <mecacell/utilities/obj3D.hpp>
#include <random>
#include <cmath>
#include "config.hpp"
#include "neurons.hpp"

template <typename cell_t, typename config_t> class Scenario;

template <typename cell_t, typename config_t> struct contractPlugin {
	Scenario<cell_t, config_t>* scenario;

	contractPlugin(Scenario<cell_t, config_t>* s) : scenario(s) {}
	template <typename W> void preBehaviorUpdate(W* w) {
		if (w->getNbUpdates() % scenario->config.t_contract == 0) {
			int ccount = 0;
			for (size_t i = 0; i < w->cells.size(); ++i) {
				if (w->cells[i]->contract) {
					ccount++;
					w->cells[i]->startContracting();
				}
			}
			MecaCell::logger<MecaCell::DBG>("Contract :: ", ccount);
		}
	}
};

template <typename cell_t, typename config_t> struct rewardPlugin {
	Scenario<cell_t, config_t>* scenario;
	MecaCell::Vec com = MecaCell::Vec(0.0, 0.0, 0.0);
	double maxd = 0.0;
	double rsignal = 0.0;
	vector<double> reward;
	vector<MecaCell::Vec> r_points;
	int r_index=0;
	std::normal_distribution<> distances;
	std::normal_distribution<> angles;

	rewardPlugin(Scenario<cell_t, config_t>* s) : scenario(s), angles(0, M_PI/8),
		distances(s->config.reward_distance_mean, s->config.reward_distance_std) {
		for (int i = 0; i < 10; i++) reward.push_back(0.0);
	}

	template <typename W> void init(W* w) {
		for (auto& c : w->cells) com += c->getPosition();
		com /= w->cells.size();
		for (auto& c : w->cells)
			if ((com-c->getPosition()).length() > maxd) maxd = (com-c->getPosition()).length();

		r_points.push_back(MecaCell::Vec(com.coords[0], com.coords[1], com.coords[2]));
	}

	template <typename W> void preBehaviorUpdate(W* w) {
		if ((w->getNbUpdates() > 0) && (w->getNbUpdates() % scenario->config.t_reward == 0)) {
			int t = w->getNbUpdates() / scenario->config.t_reward;
			MecaCell::Vec ncom(0.0, 0.0, 0.0);
			for (auto& c : w->cells) {
				ncom += c->getPosition();
				if ((com-c->getPosition()).length() > maxd) maxd = (com-c->getPosition()).length();
			}
			ncom /= w->cells.size();
			double pdistance = (r_points[r_index]-com).length()/maxd;
			double distance = (r_points[r_index]-ncom).length()/maxd;
			double r = (pdistance-distance)/0.01;
			if (r_index==0) r*=-1.0;
			r=std::min(r, 1.0); r=std::max(r, -1.0);
			reward.pop_back();
			reward.insert(reward.begin(), r);
			rsignal = (std::exp(scenario->config.reward_signal)-1)*(std::exp(r)-1);
			com = ncom;
			MecaCell::logger<MecaCell::INF>("Reward :: ", distance, ",", reward[0], ",", t, ",", r_index,
																			",", com.coords[0], ",", com.coords[1], ",", com.coords[2]);
		}
	}
};

template <typename cell_t, typename config_t> struct neuronPlugin {
	Scenario<cell_t, config_t>* scenario;
	std::uniform_real_distribution<> dis;
	SNN snn;

	neuronPlugin(Scenario<cell_t, config_t>* s) : scenario(s), dis(0, 1),
		snn(s->config.input_signal, s->config.vt, s->config.vr, s->config.aplus, s->config.aminus) {}

	template <typename W> void init(W* w) {
		MecaCell::Vec zvec(0.0, 0.0, 0.0);
		for (int i = 0; i < scenario->config.ninput; i++) {
			bool reward = dis(scenario->gen) < scenario->config.input_reward_chance;
			snn.neurons.push_back(Neuron(scenario->config.vr, NTinput, 0, reward, zvec));
		}
		for (int i = 0; i < scenario->config.nhidden; i++) {
			bool reward = dis(scenario->gen) < scenario->config.hidden_reward_chance;
			snn.neurons.push_back(Neuron(scenario->config.vr, NThidden, 0, reward, zvec));
		}
		for (size_t i = 0; i < w->cells.size(); i++) {
			snn.neurons.push_back(Neuron(scenario->config.vr, NToutput, i, false, w->cells[i]->getPosition()));
		}

		std::normal_distribution<> weight_d(scenario->config.excitatory_mean, scenario->config.excitatory_std);

		for (size_t i=0; i<snn.neurons.size(); i++) {
			vector<double> s;
			for (size_t j=0; j<snn.neurons.size(); j++) {
				if ((i != j) && ((snn.neurons[j].type == NTinput && snn.neurons[i].type == NThidden) ||
												 (snn.neurons[j].type == NThidden && snn.neurons[i].type == NThidden) ||
												 (snn.neurons[j].type == NThidden && snn.neurons[i].type == NToutput))) {
					double w = weight_d(scenario->gen);
					if (w < 0) w = 0.001; if (w > 1) w = 0.999;
					s.push_back(w);
				} else {
					s.push_back(0.0);
				}
			}
			snn.synapses.push_back(s);
		}
	}

	template <typename W> void preBehaviorUpdate(W* w) {
		if ((w->getNbUpdates() > 0) && (w->getNbUpdates() % scenario->config.t_fire == 0)) {
			int t = w->getNbUpdates() / scenario->config.t_fire;
			snn.fire(t, scenario->rp.rsignal);
			for (const auto& n : snn.neurons) {
				if (n.t_fired == t && n.type == NToutput) scenario->world.cells[n.iid]->contract = true;
			}
		}

		if ((w->getNbUpdates() > 0) && (w->getNbUpdates() % scenario->config.t_train == 0)) {
			for (auto& n : snn.neurons)
				if (n.type == NToutput) n.position = scenario->world.cells[n.iid]->getPosition();
			snn.dopamine_release(scenario->rp.reward, scenario->rp.com, scenario->rp.maxd,
													 scenario->config.dopamine_absorption, scenario->config.dopamine_delay,
													 scenario->config.dopamine_physical_attenuation);
			snn.train();
			double maxdop = 0.0;
			for (const auto& n : snn.neurons) if (n.type==NToutput && n.dopamine > maxdop) maxdop = n.dopamine;
			if (maxdop > 0.0) {
				for (const auto& n : snn.neurons) {
					if (n.type == NToutput) {
						scenario->world.cells[n.iid]->setColorHSV(n.dopamine/maxdop*300, 0.8, 0.8);
					}
				}
			}
		}
	}
};

template <typename cell_t, typename config_t> class Scenario {
 public:
	using world_t = MecaCell::World<cell_t>;
	config_t& config;
	world_t world;
	std::random_device rd;
	std::mt19937 gen;
	contractPlugin<cell_t, config_t> cp;
	rewardPlugin<cell_t, config_t> rp;
	neuronPlugin<cell_t, config_t> np;
	double duration = 0.0;

 protected:
	double currentTime = 0;
	MecaCell::SimplifiedFluidPlugin<cell_t> sfp;	// fluid dynamics

	std::vector<MecaCell::Vec> readCellCoordinates(std::string path) {
		std::vector<MecaCell::Vec> res;
		std::ifstream file(path);
		if (!file.is_open()) throw std::runtime_error("Unable to open shape file");
		std::string line;
		while (std::getline(file, line)) {
			auto vs = MecaCell::splitStr(line, ' ');
			res.push_back(MecaCell::Vec(stod(vs[0]), stod(vs[1]), stod(vs[2])));
		}
		return res;
	}

 public:
	Scenario(config_t& c) : config(c), rp(this), np(this), cp(this), gen(rd()) {}

	void init() {

		gen.seed(config.seed);
		sfp.fluidDensity = config.fluid_density;
		sfp.fluidVelocity = {0, 0, 0};
		world.registerPlugins(rp, np, cp, sfp);
		world.setDt(config.sim_dt);
		duration = config.sim_duration;

		auto cellsCoordinates = readCellCoordinates(config.sim_shape);

		for (const auto& p : cellsCoordinates) {
			cell_t* c = new cell_t(p);
			c->getBody().setRestRadius(config.cell_radius);
			c->getBody().setStiffness(config.cell_stiffness);
			c->getBody().setMass(config.cell_mass);
			c->adhCoef = config.cell_adhesion;
			world.addCell(c);
		}

		world.update();
		world.update();

		// unbreakable initial bonds & no new adhesions, only collisions
		for (auto& conn : world.cellPlugin.connections) conn.second->unbreakable = true;
		for (auto& c : world.cells) c->adhCoef = 0.0;

		rp.init(&world);
		np.init(&world);
	}

	void loop() {
		currentTime += world.getDt();
		world.update();
	}

	world_t& getWorld() { return world; }
	bool finished() { return currentTime > duration; }
};

#endif
