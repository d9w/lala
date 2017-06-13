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

template <typename cell_t, typename config_t> struct rewardPlugin {
	Scenario<cell_t, config_t>* scenario;
	MecaCell::Vec com = MecaCell::Vec(0.0, 0.0, 0.0);
	double maxd = 0.0;
	vector<double> reward;
	vector<MecaCell::Vec> r_points;
	int r_index=0;
	std::normal_distribution<> distances;
	std::normal_distribution<> angles;

	rewardPlugin(Scenario<cell_t, config_t>* s) : scenario(s), angles(0, M_PI/8),
		distances(s->config.reward_distance_mean, s->config.reward_distance_std) {
		for (int i = 0; i < 10; i++) reward.push_back(0.0);
	}

	template <typename W> void init(W* w, MecaCell::Vec com, double maxd) {
		r_points.push_back(MecaCell::Vec(com.coords[0], com.coords[1], com.coords[2]));
	}

	template <typename W> void preBehaviorUpdate(W* w) {
		if ((w->getNbUpdates() > 0) && (w->getNbUpdates() % scenario->config.t_reward == 0)) {
      // MecaCell::logger<MecaCell::DBG>("Calculating reward");
			int t = w->getNbUpdates() / scenario->config.t_reward;
			MecaCell::Vec ncom(0.0, 0.0, 0.0);
			for (auto& c : w->cells) {
				ncom += c->getPosition();
				if ((com-c->getPosition()).length() > maxd) maxd = (com-c->getPosition()).length();
			}
			ncom /= w->cells.size();
			double pdistance = (r_points[r_index]-com).length()/maxd;
			double distance = (r_points[r_index]-ncom).length()/maxd;
			double r = (pdistance-distance)*50.0;
			if (r_index==0) r*=-1.0;
			r=max(0.0, min(1.0, r));
			reward.pop_back();
			reward.insert(reward.begin(), r);
			com = ncom;
			MecaCell::logger<MecaCell::INF>("Reward :: ", distance, ",", reward[0], ",", t,
																			",", com.coords[0], ",", com.coords[1], ",", com.coords[2]);
		}
	}
};

template <typename cell_t, typename config_t> struct neuronPlugin {
	Scenario<cell_t, config_t>* scenario;
	std::uniform_real_distribution<> dis;
  std::vector<std::vector<int>> neural_map;
  int stimulus_timing = 0;
	SNN snn;

	neuronPlugin(Scenario<cell_t, config_t>* s) : scenario(s), dis(0, 1),
		snn(s->config.vt, s->config.vr, s->config.aplus, s->config.aminus) {}

	template <typename W> void init(W* w, std::vector<MecaCell::Vec> neuralCoordinates,
                                  MecaCell::Vec com, double maxd) {
		MecaCell::Vec zvec(0.0, 0.0, 0.0);
		for (int i = 0; i < scenario->config.ninput; i++) {
			snn.neurons.push_back(Neuron(scenario->config.vr, NTinput, zvec));
      std::vector<int> empty;
      neural_map.push_back(empty);
		}
		for (int i = 0; i < scenario->config.nhidden; i++) {
			snn.neurons.push_back(Neuron(scenario->config.vr, NThidden, zvec));
      std::vector<int> empty;
      neural_map.push_back(empty);
		}
		for (const auto& p : neuralCoordinates) {
      auto neuron = Neuron(scenario->config.vr, NToutput, p);
      double da = scenario->config.dopamine_absorption;
      double dd = scenario->config.dopamine_delay;
      double dpa = scenario->config.dopamine_physical_attenuation;
      double dist = (p-com).length()/maxd;
      double dtx = dd * dist * 10;
      neuron.decay = exp(-dpa * dist);
      neuron.delay = dtx;
      std::vector<int> connected_cells;
      for (size_t i = 0; i < w->cells.size(); i++) {
        double dist = (w->cells[i]->getPosition() - p).length();
        if (dist < scenario->config.neural_radius) {
          connected_cells.push_back(i);
        }
      }
      snn.neurons.push_back(neuron);
      neural_map.push_back(connected_cells);
		}

		for (size_t i=0; i<snn.neurons.size(); i++) {
			vector<double> s;
			vector<double> sd;
      if (dis(scenario->gen) < scenario->config.inhibitory_ratio) {
        snn.neurons[i].inhibitory = true;
      }
			for (size_t j=0; j<snn.neurons.size(); j++) {
        sd.push_back(0.0);
        double conn = dis(scenario->gen);
        if (conn < scenario->config.connection_ratio) {
          if (snn.neurons[i].inhibitory) {
            s.push_back(-1.0);
          } else {
            s.push_back(1.0);
          }
        } else {
          s.push_back(0.0);
        }
			}
			snn.synapses.push_back(s);
			snn.synapses_delta.push_back(sd);
		}
	}

	template <typename W> void preBehaviorUpdate(W* w) {
		if ((w->getNbUpdates() > 0) && (w->getNbUpdates() % scenario->config.t_fire == 0)) {
      // MecaCell::logger<MecaCell::DBG>("Firing");
      vector<double> inputs;
      for (int i = 0; i < (snn.neurons.size()); i++) {
        inputs.push_back(13*(dis(scenario->gen)-0.5));
      }
      if (stimulus_timing > scenario->config.stimulus_interval) {
        MecaCell::logger<MecaCell::DBG>("Stimulus");
        for (int i = 0; i < scenario->config.ninput; i++) {
          inputs[i] += scenario->config.stimulus_signal;
        }
        stimulus_timing = 0;
      }
      stimulus_timing += 1;
      // MecaCell::logger<MecaCell::DBG>("Calling fire");
			snn.fire(inputs, scenario->config.reward_signal);
      // MecaCell::logger<MecaCell::DBG>("Contraction");
      for (size_t i = 0; i < snn.neurons.size(); i++) {
				if (snn.neurons[i].fired && snn.neurons[i].type == NToutput) {
          for (size_t j = 0; j < neural_map[i].size(); j++) {
            scenario->world.cells[neural_map[i][j]]->startContracting();
          }
        }
			}
		}

		if ((w->getNbUpdates() > 0) &&
        (w->getNbUpdates() % scenario->config.t_train == 0)) {
      snn.train(scenario->config.da_factor);
    }

		if ((w->getNbUpdates() > 0) &&
        (w->getNbUpdates() % scenario->config.t_reward== 0)) {
      // MecaCell::logger<MecaCell::DBG>("Dopamine release");
			snn.dopamine_release(scenario->rp.reward,
                           scenario->config.dopamine_absorption,
                           scenario->config.dopamine_delay,
                           scenario->config.dopamine_physical_attenuation);
      double maxdop = 0.0;
      for (const auto& n : snn.neurons) {
        if (n.type==NToutput && n.dopamine > maxdop) maxdop = n.dopamine;
      }
      if (maxdop > 0.0) {
        for (size_t i = 0; i < snn.neurons.size(); i++) {
          if (snn.neurons[i].type == NToutput) {
            for (size_t j = 0; j < neural_map[i].size(); j++) {
              scenario->world.cells[neural_map[i][j]]->
                setColorHSV(snn.neurons[i].dopamine/maxdop*300, 0.8, 0.8);
            }
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
	Scenario(config_t& c) : config(c), rp(this), np(this), gen(rd()) {}

	void init() {

		gen.seed(config.seed);
		sfp.fluidDensity = config.fluid_density;
		sfp.fluidVelocity = {0, 0, 0};
		world.registerPlugins(rp, np, sfp);
		world.setDt(config.sim_dt);
		duration = config.sim_duration;

		MecaCell::logger<MecaCell::DBG>("Cell coords :: ", config.sim_shape);
		auto cellsCoordinates = readCellCoordinates(config.sim_shape);

		MecaCell::logger<MecaCell::DBG>("Cell coords :: ", config.neural_shape);
		auto neuralCoordinates = readCellCoordinates(config.neural_shape);

		for (const auto& p : cellsCoordinates) {
			cell_t* c = new cell_t(p);
			c->getBody().setRestRadius(config.cell_radius);
			c->getBody().setStiffness(config.cell_stiffness);
			c->getBody().setMass(config.cell_mass);
			c->adhCoef = config.cell_adhesion;
      c->originalRadius = config.cell_radius;
			world.addCell(c);
		}

		MecaCell::logger<MecaCell::DBG>("Added cells");
		world.update();
		world.update();

		// unbreakable initial bonds & no new adhesions, only collisions
		for (auto& conn : world.cellPlugin.connections) conn.second->unbreakable = true;
		for (auto& c : world.cells) c->adhCoef = 0.0;

		MecaCell::logger<MecaCell::DBG>("Initializing plugins");
    MecaCell::Vec com = MecaCell::Vec(0.0, 0.0, 0.0);
		for (auto& c : world.cells) com += c->getPosition();
		com /= world.cells.size();
    double maxd = 0.0;
		for (auto& c : world.cells)
			if ((com-c->getPosition()).length() > maxd) maxd = (com-c->getPosition()).length();

		rp.init(&world, com, maxd);
		np.init(&world, neuralCoordinates, com, maxd);
		MecaCell::logger<MecaCell::DBG>("Done initializing scenario");
	}

	void loop() {
		currentTime += world.getDt();
		world.update();
	}

	world_t& getWorld() { return world; }
	bool finished() { return currentTime > duration; }
};

#endif
