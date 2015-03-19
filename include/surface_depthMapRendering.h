#ifndef _SURFACE_DEPTHMAPRENDERING_PLUGIN_H_
#define _SURFACE_DEPTHMAPRENDERING_PLUGIN_H_

#include "plugin_interaction.h"

#include "camera.h"
#include "imageCoordinates.h"

#include "dialog_surface_depthMapRendering.h"
#include "Utils/Shaders/shaderSimpleColor.h"
#include "Utils/fbo.h"
#include "Utils/chrono.h"

#include "Algo/Tiling/Surface/square.h"

#include "Algo/Export/export.h"


#include <fstream>

namespace CGoGN
{

namespace SCHNApps
{

struct MapParameters
{
	MapParameters() :
		positionVBO(NULL)
	{}

	Utils::VBO* positionVBO;

	QHash<QString, Camera*> depthCameraSet;
	QHash<QString, std::vector<GLfloat>> depthImageSet;
	QHash<QString, MapHandlerGen*> projectedMapSet;
	QHash<QString, int> decompositionLevel;
};

class Surface_DepthMapRendering_Plugin : public PluginInteraction
{
	Q_OBJECT
	Q_INTERFACES(CGoGN::SCHNApps::Plugin)

public:
	Surface_DepthMapRendering_Plugin()
	{}

	~Surface_DepthMapRendering_Plugin()
	{}

private:
	virtual bool enable();
	virtual void disable();

	virtual void draw(View* view) {}
	virtual void drawMap(View* view, MapHandlerGen* map) {}

	virtual void keyPress(View* view, QKeyEvent* event) {}
	virtual void keyRelease(View* view, QKeyEvent* event) {}
	virtual void mousePress(View* view, QMouseEvent* event) {}
	virtual void mouseRelease(View* view, QMouseEvent* event) {}
	virtual void mouseMove(View* view, QMouseEvent* event) {}
	virtual void wheelEvent(View* view, QWheelEvent* event) {}

	virtual void viewLinked(View* view) {}
	virtual void viewUnlinked(View* view) {}

private slots:
	void openDepthMapRenderingDialog();
	void closeDepthMapRenderingDialog();

	void moveDownFromDialog();
	void moveUpFromDialog();

	//SCHNApps signals
	void mapAdded(MapHandlerGen* map);
	void mapRemoved(MapHandlerGen* map);

	//MapHandler signals
	void vboRemoved(Utils::VBO* vbo);

public slots: //Python calls

	void createFBO(int width, int height);

	void changePositionVBO(const QString& view, const QString& map, const QString& vbo);

	void createCameras(const QString& mapName);
	void render(const QString& mapName, const QString& directory = "/home/blettere/Projets/Results/", bool saveData = false);
	void project2DImageTo3DSpace(const QString& mapOrigin, const QString& mapGenerated);

	bool moveDownDecomposition(const QString& mapOrigin, const QString& mapGenerated);
	bool moveUpDecomposition(const QString& mapOrigin, const QString& mapGenerated);

	bool savePointCloud(const QString& mapOrigin, const QString& mapGenerated, const QString& directory = "/home/blettere/Projets/Results/");

	bool saveMergedPointCloud(const QString& mapOrigin, const std::vector<QString>& mapNames, const QString& directory = "/home/blettere/Projets/Results/");

private:
	Dialog_Surface_DepthMapRendering* m_depthMapRenderingDialog;
	QAction* m_depthMapRenderingAction;

	QHash<MapHandlerGen*, MapParameters> m_mapParameterSet;

	CGoGN::Utils::FBO* m_depthFBO;
	CGoGN::Utils::ShaderSimpleColor* m_shader;

	bool m_draw;
};

} // namespace SCHNApps

} // namespace CGoGN

#endif
