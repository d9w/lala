#define MECACELL_TERMINAL_COLORS
#define MECACELL_LOGGER_WARN_DISABLE
#define MECACELL_LOGGER_DBG_DISABLE
#include <mecacell/mecacell.h>
#include "config.hpp"

int main(int argc, char **argv) {
	Config cfg(argc, argv);
	Config::scenario_t scenario(cfg);
	auto start = chrono::high_resolution_clock::now();
	scenario.init();
	while (!scenario.finished()) scenario.loop();
	auto end = chrono::high_resolution_clock::now();
	std::cout << chrono::duration<double, milli>(end - start).count() << " ms" << endl;

	return 0;
}
