#ifndef _SURFACE_DEPTHMAPRENDERING_PLUGIN_H_
#define _SURFACE_DEPTHMAPRENDERING_PLUGIN_H_

#include "plugin_interaction.h"

#include "camera.h"

#include "dialog_surface_depthMapRendering.h"
#include "Utils/Shaders/shaderDepth.h"
#include "Utils/fbo.h"

#include "Algo/Tiling/Surface/square.h"

#include "imageCoordinates.h"

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
	virtual void drawMap(View* view, MapHandlerGen* map);

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

	void renderFromDialog();

	//SCHNApps signals
	void mapAdded(MapHandlerGen* map);
	void mapRemoved(MapHandlerGen* map);

	//MapHandler signals
	void vboRemoved(Utils::VBO* vbo);

public slots: //Python calls

	void changePositionVBO(const QString& view, const QString& map, const QString& vbo);

	void createCameras(const QString& mapName, int number = 4);
	void render(const QString& mapName, const QString& directory = "/home/blettere/Projets/Models/DepthMaps");
	void project2DImageTo3DSpace(const QString& mapOrigin, const QString& mapGenerated);

private:
	Dialog_Surface_DepthMapRendering* m_depthMapRenderingDialog;
	QAction* m_depthMapRenderingAction;

	QHash<MapHandlerGen*, MapParameters> m_mapParameterSet;

	CGoGN::Utils::ShaderDepth* m_depthShader;
	CGoGN::Utils::FBO* m_depthFBO;

	bool m_draw;
};

} // namespace SCHNApps

} // namespace CGoGN

#endif
