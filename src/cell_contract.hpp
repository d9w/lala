#ifndef CELL_HPP
#define CELL_HPP
#include <mecacell/mecacell.h>
#include <mecacell/elasticbody.hpp>
#include <mecacell/springbody.hpp>

class Cell : public MecaCell::ConnectableCell<Cell, MecaCell::SpringBody> {
	bool contracting = false;
	double contractTime = 0.0;
  double currentRadius = 30.0;

 public:
  double originalRadius = 30.0;
	double adhCoef = 25.0;
	double contractRatio = 0.9;
	double contractDuration = 0.4;
  double step = 0.0;

	using Parent = MecaCell::ConnectableCell<Cell, MecaCell::SpringBody>;
	Cell(const MecaCell::Vector3D& v) : Parent(v) {}
	double getAdhesionWith(Cell*, MecaCell::Vec) { return adhCoef; }

	void startContracting() {
		// if (!contracting) {
			// originalRadius = this->body.getRestRadius();
		contracting = true;
		contractTime = 0.0;
			// this->body.setRestRadius(originalRadius * contractRatio);
		// }
	}

	template <typename W> void updateBehavior(W& w) {
		if (contracting) {
      currentRadius = this->body.getRestRadius();
      step = (1-contractRatio) * originalRadius * 2.0 * (w.getDt() / contractDuration);
			contractTime += w.getDt();
      double cp = max(0.0, min(1.0, abs(contractTime - (contractDuration/2.0))/
                               (contractDuration/2.0)));
      double desiredRadius = (1.0-cp*(1.0-contractRatio))*originalRadius;
      std::cerr << cp << "," << desiredRadius << " ";
      if (desiredRadius < currentRadius) {
        currentRadius += step;
        this->body.setRestRadius(currentRadius);
      } else if (desiredRadius > currentRadius) {
        currentRadius -= step;
        this->body.setRestRadius(currentRadius);
      }
		}
	}
};
#endif
