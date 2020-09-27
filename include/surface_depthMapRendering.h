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
#include <QFileDialog>

void createDirectory(const QString& filename);


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
#if CGOGN_QT_DESIRED_VERSION == 5
	Q_PLUGIN_METADATA(IID "CGoGN.SCHNapps.Plugin")
#endif

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

	void createFBOFromDialog();
	void lowerResolutionFromDialog();
	void upperResolutionFromDialog();
	void saveMergedPointCloudFromDialog();
	void saveDepthMapScreenshotFromDialog();

	//SCHNApps signals
	void mapAdded(MapHandlerGen* map);
	void mapRemoved(MapHandlerGen* map);

	//MapHandler signals
	void vboRemoved(Utils::VBO* vbo);
	void selectedCellsChanged(CellSelectorGen* cs);

public slots: //Python calls

	void openDepthMapRenderingDialog();
	void closeDepthMapRenderingDialog();

	void createFBO(int width, int height);

	void changePositionVBO(const QString& view, const QString& map, const QString& vbo);

	void createCameras(const QString& camera_directory, const QString& mapName, int nb_subdiv);
	void render(const QString& mapName, const QString& directory);
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
	void saveDepthMapScreenshot(const QString& mapName,
								const QString& directory = "/home/blettere/Projets/Results/");
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

	bool readPly(
			std::vector<float>& vertices,
			std::vector<unsigned int>& connectivity,
			std::vector<uchar>& colors,
			std::string ply_filename)
	{
		std::ifstream fileStream (ply_filename, std::ios::in | std::ios::binary);

		if(fileStream.is_open())
		{
			/* Write header */
			std::string line;

			std::getline(fileStream, line);
			if(line != "ply")
			{
				std::cerr << "Incorrect input file" << std::endl;
				exit(-1);
			}

			bool is_ascii = false;

			//Read data type (ascii, binary)
			std::getline(fileStream, line);

			while(line.find("comment") != std::string::npos)
			{
				std::getline(fileStream, line);
			}

			if(line.find("ascii") != std::string::npos)
			{
				is_ascii = true;
			}

			int nbV;

			//Read number of vertices
			std::getline(fileStream, line);

			while(line.find("comment") != std::string::npos)
			{
				std::getline(fileStream, line);
			}

			if(line.find("element vertex") != std::string::npos)
			{
				std::istringstream tmp(line);
				tmp >> line >> line >> nbV;
			}
			else
			{
				std::cerr << "Incorrect information in the header" << std::endl;
				exit(-1);
			}

			//Read properties linked to vertices

			std::getline(fileStream, line);
			std::getline(fileStream, line);
			std::getline(fileStream, line);
			std::getline(fileStream, line);

			bool with_normal = false;

			if(line == "property float32 nx")
			{
				with_normal = true;

				while(line.find("property float32") != std::string::npos)
				{
					std::getline(fileStream, line);
				}
			}

			bool with_color = false;

			while(line.find("comment") != std::string::npos)
			{
				std::getline(fileStream, line);
			}

			int nb_channels = 0;

			if(line == "property uchar red")
			{
				with_color = true;

				//Jump to the next property

				while(line.find("property uchar") != std::string::npos)
				{
					++nb_channels;
					std::getline(fileStream, line);
				}
			}

			int nb_additional_attributes = 0;

			while(line.find("property float") != std::string::npos)
			{
				++nb_additional_attributes;
				std::getline(fileStream, line);
			}

			while(line.find("comment") != std::string::npos)
			{
				std::getline(fileStream, line);
			}

			//Read number of faces
			int nbF = 0;
			if(line.find("element face") != std::string::npos)
			{
				std::istringstream tmp(line);
				tmp >> line >> line >> nbF;

				//Read properties of the faces
				std::getline(fileStream, line);
				std::getline(fileStream, line);
			}

			while(line.find("comment") != std::string::npos)
			{
				std::getline(fileStream, line);
			}

			if(line != "end_header")
			{
				std::cerr << "Incorrect information in the header" << std::endl;
				exit(-1);
			}

			//End of the header

			vertices.resize(nbV*3);
			colors.resize(nbV*nb_channels, 0);
			connectivity.resize(nbF*3, 0);

			if(is_ascii)
			{

			}
			else
			{
				//Read vertices
				if(with_color && with_normal)
				{
					std::vector<float> tmp_normal(3);
					if(nb_additional_attributes > 0)
					{
						std::vector<float> tmp(nb_additional_attributes);
						for (unsigned int idx = 0; idx < nbV; ++idx)
						{
							fileStream.read(reinterpret_cast<char*>(&(vertices[3 * idx])), sizeof(float)*3);
							fileStream.read(reinterpret_cast<char*>(tmp_normal.data()), sizeof(float)*tmp_normal.size());
							fileStream.read(reinterpret_cast<char*>(&(colors[nb_channels * idx])), sizeof(unsigned char)*nb_channels);
							fileStream.read(reinterpret_cast<char*>(tmp.data()), sizeof(float)*nb_additional_attributes);
						}
					}
					else
					{
						for (unsigned int idx = 0; idx < nbV; ++idx)
						{
							fileStream.read(reinterpret_cast<char*>(&(vertices[3 * idx])), sizeof(float)*3);
							fileStream.read(reinterpret_cast<char*>(tmp_normal.data()), sizeof(float)*tmp_normal.size());
							fileStream.read(reinterpret_cast<char*>(&(colors[nb_channels * idx])), sizeof(unsigned char)*nb_channels);
						}
					}
				}
				else if(with_normal)
				{
					std::vector<float> tmp_normal(3);
					if(nb_additional_attributes > 0)
					{
						std::vector<float> tmp(nb_additional_attributes);
						for (unsigned int idx = 0; idx < nbV; ++idx)
						{
							fileStream.read(reinterpret_cast<char*>(&(vertices[3 * idx])), sizeof(float)*3);
							fileStream.read(reinterpret_cast<char*>(tmp_normal.data()), sizeof(float)*tmp_normal.size());
							fileStream.read(reinterpret_cast<char*>(&(colors[nb_channels * idx])), sizeof(unsigned char)*nb_channels);
							fileStream.read(reinterpret_cast<char*>(tmp.data()), sizeof(float)*nb_additional_attributes);
						}
					}
					else
					{
						for (unsigned int idx = 0; idx < nbV; ++idx)
						{
							fileStream.read(reinterpret_cast<char*>(&(vertices[3 * idx])), sizeof(float)*3);
							fileStream.read(reinterpret_cast<char*>(tmp_normal.data()), sizeof(float)*tmp_normal.size());
							fileStream.read(reinterpret_cast<char*>(&(colors[nb_channels * idx])), sizeof(unsigned char)*nb_channels);
						}
					}
				}
				else if(with_color)
				{
					if(nb_additional_attributes > 0)
					{
						std::vector<float> tmp(nb_additional_attributes);
						for (unsigned int idx = 0; idx < nbV; ++idx)
						{
							fileStream.read(reinterpret_cast<char*>(&(vertices[3 * idx])), sizeof(float)*3);
							fileStream.read(reinterpret_cast<char*>(&(colors[nb_channels * idx])), sizeof(unsigned char)*nb_channels);
							fileStream.read(reinterpret_cast<char*>(tmp.data()), sizeof(float)*nb_additional_attributes);
						}
					}
					else
					{
						for (unsigned int idx = 0; idx < nbV; ++idx)
						{
							fileStream.read(reinterpret_cast<char*>(&(vertices[3 * idx])), sizeof(float)*3);
							fileStream.read(reinterpret_cast<char*>(&(colors[nb_channels * idx])), sizeof(unsigned char)*nb_channels);
						}
					}
				}
				else
				{
					for (unsigned int idx = 0; idx < nbV; ++idx)
					{
						fileStream.read(reinterpret_cast<char*>(&(vertices[3 * idx])), sizeof(float)*3);
					}
				}

				//Read faces
				unsigned char face_degree;
				/* Write connectivity */
				for (unsigned int idx = 0; idx < nbF; ++idx)
				{
					fileStream.read(reinterpret_cast<char*>(&face_degree), sizeof(unsigned char));
					fileStream.read(reinterpret_cast<char*>(&(connectivity[face_degree * idx])), sizeof(int)*face_degree);
				}
			}

			vertices.shrink_to_fit();
			colors.shrink_to_fit();
			connectivity.shrink_to_fit();

			/* Close the stream */
			fileStream.close ();
		}
		return true;
	}

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
