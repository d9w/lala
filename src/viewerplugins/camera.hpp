#ifndef VIEWEREXTENSIONSCAMERA_HPP
#define VIEWEREXTENSIONSCAMERA_HPP
#include <mecacell/mecacell.h>
#include <mecacellviewer/viewer.h>

using namespace MecacellViewer;

class CenteredCameraPlugin {
 public:
	const double stiffnessConstraint = 12.0;

	const double desiredDistance = 3200.0;
	QVector3D centerOfMass = QVector3D(0, 0, 0);

	template <typename R> void onLoad(R* r) {
		r->getCamera().setMode(Camera::centered);
	}
	template <typename R> void preLoop(R* r) {
		centerOfMass = QVector3D(0, 0, 0);
		for (auto& c : r->getScenario().getWorld().cells)
			centerOfMass += toQV3D(c->getPosition());
		centerOfMass =
		    centerOfMass / static_cast<double>(r->getScenario().getWorld().cells.size());

		r->getCamera().setTarget(centerOfMass);
		double L = desiredDistance;
		double X = L - (r->getCamera().getPosition() - centerOfMass).length();

		// r->getCamera().force += QVector3D(900.0, 0.0, stiffnessConstaint * X); // slowly
		// rotate the camera arount the center of mass
	}
};
#endif
