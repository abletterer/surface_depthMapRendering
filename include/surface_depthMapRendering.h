#ifndef _SURFACE_DEPTHMAPRENDERING_PLUGIN_H_
#define _SURFACE_DEPTHMAPRENDERING_PLUGIN_H_

#include "plugin_interaction.h"

#include "dialog_surface_depthMapRendering.h"

#include "camera.h"
#include "imageCoordinates.h"

#include "Utils/Shaders/shaderSimpleColor.h"
#include "Utils/Shaders/shaderScalarFieldReal.h"
#include "Utils/fbo.h"
#include "Utils/chrono.h"

#include "Algo/Tiling/Surface/square.h"

#include "Algo/Export/export.h"
#include "Eigen/Eigen"

#include "Container/fakeAttribute.h"

#include <fstream>
#include <thread>

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
	QHash<QString, Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> > depthImageSet;
    QHash<QString, Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> > newDepthImageSet;
	QHash<QString, int> decompositionLevelSet;
	QHash<QString, MapHandlerGen*> projectedMapSet;
};

struct PointCorrespondance
{
	MapHandlerGen* map;
	Dart vertex;
};

enum Criteria
{
	DENSITY = 0,
	VISIBILITY = 1
} ;

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
	void lowerResolutionFromDialog();
	void upperResolutionFromDialog();
	void saveMergedPointCloudFromDialog();

	//SCHNApps signals
	void mapAdded(MapHandlerGen* map);
	void mapRemoved(MapHandlerGen* map);

	//MapHandler signals
	void vboRemoved(Utils::VBO* vbo);
	void selectedCellsChanged(CellSelectorGen* cs);

public slots: //Python calls

	void createFBO(int width, int height);

	void changePositionVBO(const QString& view, const QString& map, const QString& vbo);

	void createCameras(const QString& mapName, int nbMax = 12);
	void render(const QString& mapName);
	void project2DImageTo3DSpace(const QString& mapOrigin, const QString& mapGenerated);
	bool lowerResolution(const QString& mapOrigin, const QString& mapGenerated);
	bool upperResolution(const QString& mapOrigin, const QString& mapGenerated);

	/*
	 * Export de données
	 */
	bool savePointCloud(const QString& mapOrigin, const QString& mapGenerated,
						const QString& directory = "/home/blettere/Projets/Results/", const int criteria=DENSITY, const float radius=1);
	bool saveOriginalDepthMap(const QString& mapOrigin, const QString& mapGenerated,
							  const QString& directory = "/home/blettere/Projets/Results/");
	bool saveModifiedDepthMap(const QString& mapOrigin, const QString& mapGenerated,
							  const QString& directory = "/home/blettere/Projets/Results/", const int criteria=DENSITY, const float radius=1);
	bool saveMergedPointCloud(const QString& mapOrigin, const QStringList& mapNames,
							  const QString& directory = "/home/blettere/Projets/Results/", const int criteria=DENSITY, const float radius=1);
	void exportModelPly(const QString& mapName, const QString& directory = "/home/blettere/Projets/Results/");

	/*
	 * Recherche des correspondances de points
	 */
	void normalEstimation(const QString& mapOrigin, const QString& mapGenerated);
	void confidenceEstimation(const QString& mapOrigin, const QString& mapGenerated);
	void densityEstimation(const QString& mapOrigin, const QString& mapGenerated, const float radius);
	void findCorrespondingPoints(const QString& mapOrigin, const QString& mapGenerated, const int criteria=DENSITY);
    
    void updateDepthImages(const QString& mapOrigin);

	void regenerateMap(const QString& mapOrigin, const QString& mapGenerated);
	void deleteBackground(const QString& mapOrigin, const QString& mapGenerated);
	void removeUselessAttributes(const QString& mapGenerated);

	void verifyDepthMaps(const QString& mapOrigin, const QString& mapGenerated);

private:
	Dialog_Surface_DepthMapRendering* m_depthMapRenderingDialog;
	QAction* m_depthMapRenderingAction;

	QHash<MapHandlerGen*, MapParameters> m_mapParameterSet;

	CGoGN::Utils::FBO* m_fbo;
	CGoGN::Utils::ShaderSimpleColor* m_shaderSimpleColor;
	CGoGN::Utils::ShaderScalarFieldReal* m_shaderScalarFieldReal;

	bool m_draw;
	bool m_correspondance_done;

	MapHandlerGen* m_main_object;
};

} // namespace SCHNApps

} // namespace CGoGN

#endif
