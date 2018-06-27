#ifndef SCENARIO_HPP
#include <mecacell/mecacell.h>
#include <mecacell/utilities/utils.hpp>
#include <mecacell/utilities/obj3D.hpp>
#include <random>
#include <cmath>
#include <math.h>
#include "simplifiedfluidplugin.hpp"
#include "config.hpp"
#include "neurons.hpp"

template <typename cell_t, typename config_t> class Scenario;

template <typename cell_t, typename config_t> struct rewardPlugin {
	Scenario<cell_t, config_t>* scenario;
	MecaCell::Vec com = MecaCell::Vec(0.0, 0.0, 0.0);
	MecaCell::Vec init_com = MecaCell::Vec(0.0, 0.0, 0.0);
	double maxd = 0.0;
  double max_velocity = 0.0;
	vector<double> reward;
	std::normal_distribution<> distances;
	std::normal_distribution<> angles;

	rewardPlugin(Scenario<cell_t, config_t>* s) : scenario(s), angles(0, M_PI/8),
		distances(s->config.reward_distance_mean, s->config.reward_distance_std) {
		for (int i = 0; i < 10; i++) reward.push_back(0.0);
	}

	template <typename W> void init(W* w, MecaCell::Vec c, double m) {
    init_com = c; com = c; maxd = m;
	}

	template <typename W> void preBehaviorUpdate(W* w) {
    return;
		if ((w->getNbUpdates() > 0) &&
        (w->getNbUpdates() % scenario->config.t_reward == 0)) {
      // MecaCell::logger<MecaCell::DBG>("Calculating reward");
			int t = w->getNbUpdates() / scenario->config.t_reward;
			MecaCell::Vec ncom(0.0, 0.0, 0.0);
			for (auto& c : w->cells) {
				ncom += c->getPosition();
			}
			ncom /= w->cells.size();
			double pdistance = (init_com-com).length();
      double distance = (init_com-ncom).length();
			double velocity = distance-pdistance;
      double r = 0.0;
      if (t > 10 && velocity > max_velocity) {
        if (max_velocity > 0.0) {
          r = 10.0 * (velocity - max_velocity) / max_velocity;
        }
        max_velocity = velocity;
      }
      max_velocity *= 0.99;
			reward.pop_back();
			reward.insert(reward.begin(), r);
			com = ncom;
      scenario->np.stimulus_timing = 0;

      /* LOG INFO */
      double weight_mean = scenario->np.snn.get_weight_mean();
      double weight_std = scenario->np.snn.get_weight_std(weight_mean);
      double da_mean = scenario->np.snn.get_da_mean();
      double da_std = scenario->np.snn.get_da_std(da_mean);
			MecaCell::logger<MecaCell::INF>("Reward :: ", distance, ",", velocity,",",
                                      max_velocity, ",",
                                      reward[0], ",", t, ",", weight_mean, ",",
                                      weight_std, ",", da_mean, ",", da_std);
      /* END LOG INFO */
		}
	}
};

template <typename cell_t, typename config_t> struct neuronPlugin {
	Scenario<cell_t, config_t>* scenario;
	std::uniform_real_distribution<> dis;
  std::vector<std::vector<int>> neural_map;
  int stimulus_timing = 0;
  int stimulus_group = 0;
  double max_dopamine = 0.0;
	SNN snn;

	neuronPlugin(Scenario<cell_t, config_t>* s) : scenario(s), dis(0, 1),
		snn(s->config.vt, s->config.vr, s->config.aplus, s->config.aminus) {}

	template <typename W> void init(W* w, std::vector<MecaCell::Vec> neuralCoordinates,
                                  MecaCell::Vec com, double maxd) {
    return;
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
      double dist = (p-com).length()/maxd;
      neuron.decay = exp(-scenario->config.dopamine_physical_attenuation*dist);
      neuron.delay = scenario->config.dopamine_delay * dist * 10;
      snn.neurons.push_back(neuron);
      std::vector<int> empty;
      neural_map.push_back(empty);
		}

    for (size_t i = 0; i < w->cells.size(); i++) {
      double min_dist = 1e9;
      int min_ind = 0;
      for (size_t j=0; j<snn.neurons.size(); j++) {
        if (snn.neurons[j].type == NToutput) {
          double dist = (w->cells[i]->getPosition() -
                         snn.neurons[j].position).length();
          if (dist < min_dist) {
            min_dist = dist;
            min_ind = j;
          }
        }
      }
      neural_map[min_ind].push_back(i);
    }

		for (size_t i=0; i<snn.neurons.size(); i++) {
			vector<double> s;
			vector<double> sd;
      if (dis(scenario->gen) < scenario->config.inhibitory_ratio) {
        snn.neurons[i].inhibitory = true;
      }
			for (size_t j=0; j<snn.neurons.size(); j++) {
        sd.push_back(0.0);
        // if ((snn.neurons[i].type == NTinput &&  snn.neurons[j].type == NThidden) ||
        //     (snn.neurons[i].type == NThidden && snn.neurons[j].type == NThidden) ||
        //     (snn.neurons[i].type == NThidden && snn.neurons[j].type == NToutput)) {
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
      // } else {
      //     s.push_back(0.0);
      //   }
			}
			snn.synapses.push_back(s);
			snn.synapses_delta.push_back(sd);
		}
	}

	template <typename W> void preBehaviorUpdate(W* w) {
    return;
		if ((w->getNbUpdates() > 0) &&
        (w->getNbUpdates() % scenario->config.t_fire == 0)) {
      // MecaCell::logger<MecaCell::DBG>("Firing");
      vector<double> inputs;
      for (int i = 0; i < (snn.neurons.size()); i++) {
        inputs.push_back(13*(dis(scenario->gen)-0.5));
      }
      if ((stimulus_timing > scenario->config.stimulus_interval) &&
          dis(scenario->gen) > 0.75) {
        MecaCell::logger<MecaCell::DBG>("Stimulus");
        int groupsize = int(scenario->config.ninput /
                            scenario->config.stimulus_groups);
        int start = stimulus_group * groupsize;
        int end = (stimulus_group + 1)* groupsize;
        for (int i = start; i < end; i++) {
          inputs[i] += scenario->config.stimulus_signal;
        }
        stimulus_group = (stimulus_group + 1) %
          (scenario->config.stimulus_groups - 1);
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
      double maxda = 0.0;
      double minda = 100.0;
      for (const auto& n : snn.neurons) {
        if (n.type==NToutput && n.dopamine > maxda) maxda = n.dopamine;
        if (n.type==NToutput && n.dopamine < minda) minda = n.dopamine;
      }
      if (maxda > 0.0) {
        if (maxda > max_dopamine) max_dopamine = maxda;
        double denom = maxda - minda;
        if (denom == 0.0) denom = 1.0;
        for (size_t i = 0; i < snn.neurons.size(); i++) {
          if (snn.neurons[i].type == NToutput) {
            for (size_t j = 0; j < neural_map[i].size(); j++) {
              scenario->world.cells[neural_map[i][j]]->
                setColorHSV((snn.neurons[i].dopamine-minda)/denom*300, 0.8,
                            0.8*maxda/max_dopamine);
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
			c->getBody().setRadius(config.cell_radius);
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
