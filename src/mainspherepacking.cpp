#include <mecacell/mecacell.h>
#include <iostream>
#include <mecacell/utilities/obj3D.hpp>
#include <vector>

std::vector<MecaCell::Vec> spherePacking(double radius, double lngth, double wdth,
                                         double hght, const MecaCell::Scene3D& scene) {
	std::vector<MecaCell::Vec> res;
	for (double i = 0; i < lngth; ++i) {
		for (double j = 0; j < wdth; ++j) {
			for (double k = 0; k < hght; ++k) {
				double X = i * radius;
				double Y = j * radius;
				double Z = k * radius;
				if ((static_cast<int>(j) % 2) == 0) X = X + (radius / 2.0);
				if ((static_cast<int>(k) % 2) == 0) Y = Y + (radius / 2.0);
				MecaCell::Vec pos(X - lngth * radius * 0.5, Y - wdth * radius * 0.5,
				                  Z - hght * radius * 0.5);
				auto gr = scene.isInside(pos);
				if (gr.size() > 0) {
					res.push_back(pos);
					// MecaCell::logger<MecaCell::DBG>(pos, " is inside");
				} else {
					// MecaCell::logger<MecaCell::DBG>(pos, " is outside, gr = ", gr.size());
				}
			}
		}
	}
	return res;
}

int main(int argc, char** argv) {
	if (argc != 3) {
		std::cerr << " usage: " << argv[0] << " path/to/shape.obj distBetweenCells"
		          << std::endl;
		return 1;
	}
	MecaCell::Scene3D shape(argv[1]);
	shape.scale(MecaCell::Vec(50, 50, 50));
	double distBetweenCells = stod(argv[2]);
	auto cellsCoordinates = spherePacking(distBetweenCells, 100, 100, 100, shape);
	for (auto& p : cellsCoordinates)
		std::cerr << p.x() << " " << p.y() << " " << p.z() << std::endl;
	return 0;
}
