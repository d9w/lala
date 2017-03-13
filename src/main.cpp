#include "include/mecacell/mecacell/mecacell.h"
#include "cell.h"
#include "scenario.h"

// int main(int argc, char **argv) {
//   MecacellViewer::Viewer<MyScenario> v;
//   return v.exec(argc, argv);
// }

int main(int argc, char **argv) {
  MyScenario scenario;
  for (auto step=0; step<10; step++) {
    scenario.loop();
  }
  return 0;
}
