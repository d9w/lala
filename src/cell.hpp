#ifndef CELL_HPP
#define CELL_HPP
#include <mecacell/mecacell.h>
#include <mecacell/elasticbody.hpp>
#include <mecacell/springbody.hpp>

class Cell : public MecaCell::ConnectableCell<Cell, MecaCell::SpringBody> {
	bool contracting = false;
	double contractTime = 0.0;

 public:
	double originalRadius = 30.0;
	double adhCoef = 25.0;
	double contractRatio = 0.9;
	double contractDuration = 0.4;

	using Parent = MecaCell::ConnectableCell<Cell, MecaCell::SpringBody>;
	Cell(const MecaCell::Vector3D& v) : Parent(v) {}
	double getAdhesionWith(Cell*, MecaCell::Vec) { return adhCoef; }

	void startContracting() {
		contractTime = 0.0;
    if (!contracting) {
      contracting = true;
			this->body.setRestRadius(originalRadius * contractRatio);
		}
	}

	template <typename W> void updateBehavior(W& w) {
		if (contracting) {
			contractTime += w.getDt();
			if (contractTime > contractDuration) {
				contracting = false;
				contractTime = 0.0;
				this->body.setRestRadius(originalRadius);
			}
		}
	}
};
#endif
