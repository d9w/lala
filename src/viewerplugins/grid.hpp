#ifndef DRAWEXTENSIONS_HPP
#define DRAWEXTENSIONS_HPP
#include <mecacellviewer/viewer.h>
#include <mecacellviewer/primitives/sphere.hpp>
#include <mecacell/utilities/grid.hpp>

struct GridViewer {
	QOpenGLShaderProgram shader;
	Cube cube;

 public:
	void load() {
		shader.addShaderFromSourceCode(QOpenGLShader::Vertex,
		                               shaderWithHeader(":/shaders/mvp.vert"));
		shader.addShaderFromSourceCode(QOpenGLShader::Fragment,
		                               shaderWithHeader(":/shaders/flat.frag"));
		shader.link();
		cube.load(shader);
	}

	template <typename R> void call(R *r) {
		MecaCell::Grid<typename R::Cell *> g(38);
		for (auto &c : r->getScenario().getWorld().cells) g.insert(c);

		const QMatrix4x4 &view = r->getViewMatrix();
		const QMatrix4x4 &projection = r->getProjectionMatrix();

		shader.bind();
		cube.vao.bind();
		shader.setUniformValue(shader.uniformLocation("projection"), projection);
		shader.setUniformValue(shader.uniformLocation("view"), view);

		double cellSize = g.getCellSize();

		for (const auto &c : g.getContent()) {
			QMatrix4x4 model;
			double v = g.vecToColor(c.first);
			double h = v / 8.0;
			QColor col = QColor::fromHsvF(h, 1, 0.65);
			model.translate(toQV3D(c.first) * cellSize + 0.5 * cellSize * QVector3D(1, 1, 1));
			model.scale(cellSize * 0.5, cellSize * 0.5, cellSize * 0.5);
			QMatrix4x4 nmatrix = (model).inverted().transposed();
			shader.setUniformValue(shader.uniformLocation("model"), model);
			shader.setUniformValue(shader.uniformLocation("normalMatrix"), nmatrix);
			shader.setUniformValue(shader.uniformLocation("color"), col);
			GL->glDrawElements(GL_TRIANGLES, cube.indices.size(), GL_UNSIGNED_INT, 0);
		}
		cube.vao.release();
		shader.release();
	}
};

class GridViewerPlugin {
	GridViewer gv;

 public:
	template <typename R> void onLoad(R *renderer) {
		MenuElement<R> *nativeDisplayMenu = renderer->getDisplayMenu();
		MenuElement<R> gridviewer = {"Grid Viewer", false};
		gv.load();
		gridviewer.onToggled = [&](R *r, MenuElement<R> *me) {
			if (me->isChecked())
				r->addPaintStepsMethods(42, [&](R *r2) { gv.call(r2); });
			else
				r->erasePaintStepsMethods(42);
		};
		nativeDisplayMenu->at("Cells").add(gridviewer);
	}
};

struct PointViewer {
	QOpenGLShaderProgram shader;
	IcoSphere sphere;

 public:
	PointViewer() : sphere(1) {}

	void load() {
		shader.addShaderFromSourceCode(QOpenGLShader::Vertex,
		                               shaderWithHeader(":/shaders/mvp.vert"));
		shader.addShaderFromSourceCode(QOpenGLShader::Fragment,
		                               shaderWithHeader(":/shaders/flat.frag"));
		shader.link();
		sphere.load(shader);
	}

	template <typename R> void call(R *r) {
		const QMatrix4x4 &view = r->getViewMatrix();
		const QMatrix4x4 &projection = r->getProjectionMatrix();

		shader.bind();
		sphere.vao.bind();
		shader.setUniformValue(shader.uniformLocation("projection"), projection);
		shader.setUniformValue(shader.uniformLocation("view"), view);

		double radius = r->getScenario().rp.maxd;
		for (const auto &p : r->getScenario().rp.r_points) {
			QMatrix4x4 model;
			QColor col = QColor::fromHsvF(0.6, 1, 0.95);
			col.setAlphaF(0.16);
			model.translate(toQV3D(p));
			model.scale(radius);
			QMatrix4x4 nmatrix = (model).inverted().transposed();
			shader.setUniformValue(shader.uniformLocation("model"), model);
			shader.setUniformValue(shader.uniformLocation("normalMatrix"), nmatrix);
			shader.setUniformValue(shader.uniformLocation("color"), col);
			GL->glDrawElements(GL_TRIANGLES, sphere.indices.size(), GL_UNSIGNED_INT, 0);
		}
		sphere.vao.release();
		shader.release();
	}
};

class PointViewerPlugin {
	PointViewer gv;

 public:
	template <typename R> void onLoad(R *renderer) {
		MenuElement<R> *nativeDisplayMenu = renderer->getDisplayMenu();
		MenuElement<R> gridviewer = {"Points Viewer", true};
		gv.load();
		gridviewer.onToggled = [&](R *r, MenuElement<R> *me) {
			if (me->isChecked())
				r->addPaintStepsMethods(60, [&](R *r2) { gv.call(r2); });
			else
				r->erasePaintStepsMethods(60);
		};
		nativeDisplayMenu->at("Cells").add(gridviewer);
	}
};
#endif
