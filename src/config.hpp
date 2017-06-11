#ifndef LALACONFIG_HPP
#define LALACONFIG_HPP
#include <mecacell/mecacell.h>
#include <fstream>
#include "cell.hpp"
#include "external/cxxopts/cxxopts.hpp"
#include "external/json/json.hpp"
#include "scenario.hpp"

#define CHKPARAM(paramName)         \
	if (it.key() == "" #paramName "") \
	paramName = it.value(),           \
	MecaCell::logger<MecaCell::INF>("Config :: \"", it.key(), "\" = ", it.value())

template <typename Cell, typename Conf> class Scenario;
struct Config {
	using json = nlohmann::json;

	// ---------	 STATIC CONFIG	----------
	using scenario_t = Scenario<Cell, Config>;

	// --------		DYNAMIC CONFIG	----------
	// params and their default values
	double sim_duration = 1000.0;
	double sim_dt = 0.1;
	int t_fire = 100;
	int t_train = 500;
	int t_contract = 500;
	int t_reward = 500;
	std::string sim_shape = "";
	std::string neural_shape = "";
	double cell_radius = 10.0;
	double cell_mass = 0.001;
	double cell_stiffness = 10.0;
	double cell_adhesion = 5.0;
	double aminus = -0.00375;
	double aplus = 0.005;
	double vt = 64.0;
	double vr = 0.0;
	double excitatory_mean = 0.8;
	double excitatory_std = 0.05;
	double input_signal = 1.0;
	double reward_signal = 10.0;
	int ninput = 5;
	int nhidden = 50;
  double neural_radius = 200;
	double input_reward_chance = 0.2;
	double hidden_reward_chance = 0.1;
	double dopamine_absorption = 0.1;
	double dopamine_physical_attenuation = 0.0;
	double dopamine_delay = 0.0;
	double fluid_density = 1e-4;
	double reward_distance_mean = 5.0;
	double reward_distance_std = 2.0;
	int seed = 0;

	Config(int argc, char** argv) {
		cxxopts::Options options("lala", "lala experiment program");
		options.add_options()("f,file", "configuration file", cxxopts::value<std::string>());
		options.parse(argc, argv);

		if (!options.count("file"))
			MecaCell::logger<MecaCell::WARN>("No configuration file specified. Using defaults");
		else
			load(options["file"].as<std::string>());
	}

	// loads a conf file
	void load(std::string file) {
		std::ifstream t(file);
		std::stringstream buffer;
		buffer << t.rdbuf();
		auto o = json::parse(buffer.str());
		for (auto it = o.begin(); it != o.end(); ++it) {
			CHKPARAM(sim_duration);
			else CHKPARAM(sim_dt);
			else CHKPARAM(t_fire);
			else CHKPARAM(t_train);
			else CHKPARAM(t_contract);
			else CHKPARAM(t_reward);
			else CHKPARAM(sim_shape);
			else CHKPARAM(neural_shape);
			else CHKPARAM(cell_radius);
			else CHKPARAM(cell_mass);
			else CHKPARAM(cell_stiffness);
			else CHKPARAM(cell_adhesion);
			else CHKPARAM(aminus);
			else CHKPARAM(aplus);
			else CHKPARAM(vt);
			else CHKPARAM(vr);
			else CHKPARAM(excitatory_mean);
			else CHKPARAM(excitatory_std);
			else CHKPARAM(input_signal);
			else CHKPARAM(reward_signal);
			else CHKPARAM(ninput);
			else CHKPARAM(nhidden);
			else CHKPARAM(input_reward_chance);
			else CHKPARAM(hidden_reward_chance);
			else CHKPARAM(neural_radius);
			else CHKPARAM(dopamine_absorption);
			else CHKPARAM(dopamine_physical_attenuation);
			else CHKPARAM(dopamine_delay);
			else CHKPARAM(fluid_density);
			else CHKPARAM(reward_distance_mean);
			else CHKPARAM(reward_distance_std);
			else CHKPARAM(seed);
			else MecaCell::logger<MecaCell::WARN>("Config :: unknown field \"", it.key(), "\"");
		}
	}
};
#endif
