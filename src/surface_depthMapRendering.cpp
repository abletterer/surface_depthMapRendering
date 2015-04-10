#include "surface_depthMapRendering.h"

#include "mapHandler.h"

namespace CGoGN
{

namespace SCHNApps
{

bool Surface_DepthMapRendering_Plugin::enable()
{
	m_depthMapRenderingDialog = new Dialog_Surface_DepthMapRendering(m_schnapps);

	m_depthMapRenderingAction = new QAction("Depth-map rendering", this);

	m_schnapps->addMenuAction(this, "Surface;Depth-map rendering", m_depthMapRenderingAction);

	m_draw = false;
	m_correspondance_done = false;

	m_shaderSimpleColor = new CGoGN::Utils::ShaderSimpleColor();
	m_shaderScalarFieldReal = new CGoGN::Utils::ShaderScalarFieldReal();

	registerShader(m_shaderSimpleColor);
	registerShader(m_shaderScalarFieldReal);

	connect(m_depthMapRenderingAction, SIGNAL(triggered()), this, SLOT(openDepthMapRenderingDialog()));

	connect(m_schnapps, SIGNAL(mapAdded(MapHandlerGen*)), this, SLOT(mapAdded(MapHandlerGen*)));
	connect(m_schnapps, SIGNAL(mapRemoved(MapHandlerGen*)), this, SLOT(mapRemoved(MapHandlerGen*)));

	foreach(MapHandlerGen* map, m_schnapps->getMapSet().values())
		mapAdded(map);

	return true;
}

void Surface_DepthMapRendering_Plugin::disable()
{
	if(m_shaderSimpleColor)
	{
		delete m_shaderSimpleColor;
	}
	if(m_shaderScalarFieldReal)
	{
		delete m_shaderScalarFieldReal;
	}
	if(m_depthFBO)
	{
		delete m_depthFBO;
	}

	disconnect(m_depthMapRenderingAction, SIGNAL(triggered()), this, SLOT(openDepthMapRenderingDialog()));

	disconnect(m_schnapps, SIGNAL(mapAdded(MapHandlerGen*)), this, SLOT(mapAdded(MapHandlerGen*)));
	disconnect(m_schnapps, SIGNAL(mapRemoved(MapHandlerGen*)), this, SLOT(mapRemoved(MapHandlerGen*)));
}

void Surface_DepthMapRendering_Plugin::openDepthMapRenderingDialog()
{
	m_depthMapRenderingDialog->show();
}

void Surface_DepthMapRendering_Plugin::closeDepthMapRenderingDialog()
{
	m_depthMapRenderingDialog->close();
}

void Surface_DepthMapRendering_Plugin::moveDownFromDialog()
{
	QList<QListWidgetItem*> currentItems = m_depthMapRenderingDialog->list_maps->selectedItems();
	QList<QListWidgetItem*> currentItemsGenerated = m_depthMapRenderingDialog->list_maps_generated->selectedItems();
	if(!currentItems.empty() && !currentItemsGenerated.empty())
	{
		moveDownDecomposition(currentItems[0]->text(), currentItemsGenerated[0]->text());
	}
}

void Surface_DepthMapRendering_Plugin::moveUpFromDialog()
{
	QList<QListWidgetItem*> currentItems = m_depthMapRenderingDialog->list_maps->selectedItems();
	QList<QListWidgetItem*> currentItemsGenerated = m_depthMapRenderingDialog->list_maps_generated->selectedItems();
	if(!currentItems.empty() && !currentItemsGenerated.empty())
	{
		moveUpDecomposition(currentItems[0]->text(), currentItemsGenerated[0]->text());
	}
}

void Surface_DepthMapRendering_Plugin::mapAdded(MapHandlerGen *map)
{
	connect(map, SIGNAL(vboRemoved(Utils::VBO*)), this, SLOT(vboRemoved(Utils::VBO*)));
}

void Surface_DepthMapRendering_Plugin::mapRemoved(MapHandlerGen *map)
{
	disconnect(map, SIGNAL(vboRemoved(Utils::VBO*)), this, SLOT(vboRemoved(Utils::VBO*)));
}

void Surface_DepthMapRendering_Plugin::vboRemoved(Utils::VBO *vbo)
{
	QHash<MapHandlerGen*, MapParameters>::iterator i;
	for (i = m_mapParameterSet.begin(); i != m_mapParameterSet.end(); ++i)
	{
		MapParameters& mapParams = i.value();
		if(mapParams.positionVBO == vbo)
		{
			mapParams.positionVBO = NULL;
		}
	}
}

void Surface_DepthMapRendering_Plugin::createFBO(int width, int height)
{
	m_depthFBO = new CGoGN::Utils::FBO(width, height);
	m_depthFBO->createAttachDepthTexture();
	m_depthFBO->createAttachColorTexture(GL_RGBA);
}

void Surface_DepthMapRendering_Plugin::changePositionVBO(const QString& view, const QString& map, const QString& vbo)
{
	View* v = m_schnapps->getView(view);
	MapHandlerGen* m = m_schnapps->getMap(map);
	if(v && m)
	{
		Utils::VBO* vbuf = m->getVBO(vbo);
		m_mapParameterSet[m].positionVBO = vbuf;
		if(v->isSelectedView())
		{
			if(v->isLinkedToMap(m))	v->updateGL();
		}
	}
}

void Surface_DepthMapRendering_Plugin::createCameras(const QString& mapName, int nbMax)
{
	MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
	MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);

	if(mh_map && m_mapParameterSet.contains(mhg_map))
	{
		QString baseName("DepthCamera-");

		MapParameters& mapParams = m_mapParameterSet[mhg_map];

		//Vertices coordinates of icosahedron -> regular sampling of a sphere
		std::vector<qglviewer::Vec> positions;
		positions.reserve(12);
		positions.push_back(qglviewer::Vec(0,1,2));
		positions.push_back(qglviewer::Vec(0,1,-2));
		positions.push_back(qglviewer::Vec(0,-1,2));
		positions.push_back(qglviewer::Vec(0,-1,-2));

		positions.push_back(qglviewer::Vec(1,2,0));
		positions.push_back(qglviewer::Vec(1,-2,0));
		positions.push_back(qglviewer::Vec(-1,2,0));
		positions.push_back(qglviewer::Vec(-1,-2,0));

		positions.push_back(qglviewer::Vec(2,0,1));
		positions.push_back(qglviewer::Vec(2,0,-1));
		positions.push_back(qglviewer::Vec(-2,0,1));
		positions.push_back(qglviewer::Vec(-2,0,-1));

		qglviewer::Vec bb_min = mh_map->getBBmin();
		qglviewer::Vec bb_max = mh_map->getBBmax();

		qglviewer::Vec center = (bb_min+bb_max)/2.f;

		for(int i = 0; i < nbMax; ++i)
		{
			QString cameraName(baseName);
			cameraName.append(QString::number(i));
			Camera* camera = m_schnapps->addCamera(cameraName);

			qglviewer::Vec camera_position(camera->position());

			float radius = qAbs(camera_position.z - center.z);
			++radius;	//To avoid problems when camera is placed at the center of the scene

			camera_position.x = center.x + radius*positions[i].x;
			camera_position.y = center.y + radius*positions[i].y;
			camera_position.z = center.z + radius*positions[i].z;

			camera->setPosition(camera_position);

			camera->lookAt(center);

			camera->setSceneBoundingBox(bb_min,bb_max);
			camera->showEntireScene();

			QString generatedName(mapName);
			generatedName += "-" + cameraName;

			mapParams.depthCameraSet[generatedName] = camera;
		}
	}
}

void Surface_DepthMapRendering_Plugin::render(const QString& mapName)
{
	MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
	MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);

	if(m_depthFBO && mh_map && m_mapParameterSet.contains(mhg_map))
	{
		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
		if(!position.isValid())
		{
			CGoGNerr << "position attribute is not valid" << CGoGNendl;
			return;
		}

		MapParameters& mapParams = m_mapParameterSet[mhg_map];

		int width = m_depthFBO->getWidth(), height = m_depthFBO->getHeight();
		m_shaderSimpleColor->setAttributePosition(mapParams.positionVBO);

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> pixels;
		pixels.setZero(width, height);

		Utils::Chrono chrono;
		chrono.start();

		for(QHash<QString, Camera*>::iterator it = mapParams.depthCameraSet.begin(); it != mapParams.depthCameraSet.end(); ++it)
		{
			Camera* camera = it.value();

			QString cameraName(camera->getName());

			QString generatedName(mapName);
			generatedName += "-" + cameraName;

			m_schnapps->getSelectedView()->setCurrentCamera(camera, false);

			camera->setZNear(camera->zNear());
			camera->setZFar(camera->zFar());
			camera->setStandard(false);

			m_depthFBO->bind();
			glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );	//To clean the color and depth textures

			mh_map->draw(m_shaderSimpleColor, CGoGN::Algo::Render::GL2::TRIANGLES);	//Render the map into the FrameBufferObject

			//Read pixels and store them in an array
			glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, pixels.data());
			glBindTexture(GL_TEXTURE_2D, 0);
			m_depthFBO->unbind();

			m_schnapps->getSelectedView()->setCurrentCamera("camera_0", false);

			MapHandlerGen* mhg_generated = m_schnapps->addMap(generatedName, 2);
			mapParams.projectedMapSet[generatedName] = mhg_generated;

			MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);
			PFP2::MAP* generated_map = mh_generated->getMap();

			VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinatesGenerated = mh_generated->addAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
			VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinatesGenerated = mh_generated->addAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");

			Algo::Surface::Tilings::Square::Grid<PFP2> grid(*generated_map, width-1, height-1);
			grid.embedIntoGrid(planeCoordinatesGenerated, width-1, height-1);

			std::vector<Dart>& vDarts = grid.getVertexDarts();

			for(int i = 0; i < width; ++i)
			{
				for(int j = 0; j < height; ++j)
				{
					//Set plane coordinates in [-1;1]
					planeCoordinatesGenerated[vDarts[j*width+i]][0] /= (width-1)/2.f;
					planeCoordinatesGenerated[vDarts[j*width+i]][1] /= (height-1)/2.f;

					imageCoordinatesGenerated[vDarts[j*width+i]].setCoordinates(i, j);
				}
			}

			pixels.array() = pixels.array()*2-1;	//Put depth values in the range [-1;1]

			mapParams.depthImageSet[generatedName] = pixels;

			mh_generated->notifyConnectivityModification(false);
			mh_generated->notifyAttributeModification(planeCoordinatesGenerated, false);
			mh_generated->notifyAttributeModification(imageCoordinatesGenerated, false);
			project2DImageTo3DSpace(mapName, generatedName);
		}

		CGoGNout << "Temps d'échantillonnage : " << chrono.elapsed() << " ms " << CGoGNflush;
		CGoGNout << "pour " << mapParams.depthCameraSet.size() << " vue(s) différente(s) " << CGoGNflush;
		CGoGNout << "de taille " << width << "x" << height << CGoGNflush;
		CGoGNout << " sur un objet composé de " << mh_map->getMap()->getNbCells(VERTEX) << " point(s)" << CGoGNendl;

		m_schnapps->getSelectedView()->updateGL();
	}
}

void Surface_DepthMapRendering_Plugin::project2DImageTo3DSpace(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);

	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(m_depthFBO && mh_origin && mh_generated && m_mapParameterSet.contains(mhg_origin))
	{
		MapParameters& mapParams = m_mapParameterSet[mhg_origin];
		PFP2::MAP* generated_map = mh_generated->getMap();

		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("position");
		if(!position.isValid())
		{
			position = mh_generated->addAttribute<PFP2::VEC3, VERTEX>("position");
		}

		VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinates = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
		if(!planeCoordinates.isValid())
		{
			CGoGNerr << "PlaneCoordinates attribute is not valid" << CGoGNendl;
			return;
		}

		VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_generated->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
		if(!imageCoordinates.isValid())
		{
			CGoGNerr << "ImageCoordinates attribute is not valid" << CGoGNendl;
			return;
		}

		Camera* camera = mapParams.depthCameraSet[mapGenerated];

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> pixels = mapParams.depthImageSet[mapGenerated];

		GLdouble mvp_matrix[16];
		camera->getModelViewProjectionMatrix(mvp_matrix);

		PFP2::MATRIX44 model_view_projection_matrix, model_view_projection_matrix_inv;

		for(int i = 0; i < 4; ++i)
		{
			for(int j = 0; j < 4; ++j)
			{
				model_view_projection_matrix(i,j) = mvp_matrix[i+4*j];
			}
		}

		model_view_projection_matrix.invert(model_view_projection_matrix_inv);

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			float color = pixels(imageCoordinates[d].getXCoordinate(),imageCoordinates[d].getYCoordinate());

			if(fabs(1-color) > FLT_EPSILON)
			{
				PFP2::VEC4 pos = PFP2::VEC4(planeCoordinates[d][0], planeCoordinates[d][1], color, 1.f);

				pos = model_view_projection_matrix_inv*pos;

				position[d] = PFP2::VEC3(pos[0]/pos[3], pos[1]/pos[3], pos[2]/pos[3]);
			}
		}

		generated_map->removeAttribute(planeCoordinates);

		mh_generated->notifyConnectivityModification(false);
		mh_generated->notifyAttributeModification(position, false);
		mh_generated->updateBB(position);
	}
}

bool Surface_DepthMapRendering_Plugin::savePointCloud(const QString& mapOrigin, const QString& mapGenerated, const QString& directory)
{
	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated );

	if(!directory.isEmpty() && mh_generated)
	{
		PFP2::MAP* generated_map = mh_generated->getMap();

		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("position");
		if(!position.isValid())
		{
			CGoGNerr << "position attribute is not valid" << CGoGNendl;
			return false;
		}

		QString filename(directory);
		filename += "/" + mapOrigin + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += "PointClouds/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += QString::number(m_depthFBO->getWidth()) + "x" + QString::number(m_depthFBO->getHeight()) + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		if(m_correspondance_done)
		{
			filename += mapGenerated + "-Without.ply";
		}
		else
		{
			filename += mapGenerated + "-With.ply";
		}

		return Algo::Surface::Export::exportPLYVert<PFP2>(*generated_map, position, filename.toStdString().c_str(), false);
	}

	return false;
}

bool Surface_DepthMapRendering_Plugin::saveOriginalDepthMap(const QString& mapOrigin, const QString& mapGenerated, const QString& directory)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(!directory.isEmpty() && mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin))
	{
		MapParameters& mapParams = m_mapParameterSet[mh_origin];
		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> pixels = mapParams.depthImageSet[mapGenerated];
		Camera* camera = mapParams.depthCameraSet[mapGenerated];

		QString filename(directory);
		filename += "/" + mapOrigin + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += "DepthMaps/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += QString::number(m_depthFBO->getWidth()) + "x" + QString::number(m_depthFBO->getHeight()) + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += mapGenerated;

		std::ofstream out;
		out.open(filename.toStdString() + "-originalDepthMap.dat", std::ios::out);
		if(!out.good())
		{
			CGoGNerr << "Unable to open file" << CGoGNendl;
			return false;
		}

		for(int j = m_depthFBO->getHeight()-1; j >= 0; --j)
		{
			for(unsigned int i = 0; i < m_depthFBO->getWidth(); ++i)
			{
				out << pixels(i, j) << " " << std::flush;
			}
			out << std::endl;
		}

		out.close();

		out.open(filename.toStdString() + "-MVPMatrix.dat", std::ios::out);
		if(!out.good())
		{
			CGoGNerr << "Unable to open file" << CGoGNendl;
			return false;
		}

		GLdouble mvp_matrix[16];

		camera->getModelViewProjectionMatrix(mvp_matrix);

		for(int i = 0; i < 4; ++i)
		{
			for(int j = 0; j < 4; ++j)
			{
				out << mvp_matrix[i+4*j] << " " << std::flush;
			}
			out << std::endl;
		}

		out.close();

		return true;
	}

	return false;
}

bool Surface_DepthMapRendering_Plugin::saveModifiedDepthMap(const QString& mapOrigin, const QString& mapGenerated, const QString& directory)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(!directory.isEmpty() && mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin))
	{
		PFP2::MAP* generated_map = mh_generated->getMap();
		VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_generated->getAttribute<ImageCoordinates,VERTEX>("ImageCoordinates");
		if(!imageCoordinates.isValid())
		{
			CGoGNerr << "ImageCoordinates attribute is not valid" << CGoGNendl;
			return false;
		}

		int width = m_depthFBO->getWidth(), height = m_depthFBO->getHeight();

		MapParameters& mapParams = m_mapParameterSet[mh_origin];
		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> pixels = mapParams.depthImageSet[mapGenerated];
		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> pixelsLeft;
		pixelsLeft.setOnes(width, height);

		QString filename(directory);
		filename += "/" + mapOrigin + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += "DepthMaps/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += QString::number(width) + "x" + QString::number(height) + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += mapGenerated;

		std::ofstream out;
		out.open(filename.toStdString() + "-modifiedDepthMap.dat", std::ios::out);
		if(!out.good())
		{
			CGoGNerr << "Unable to open file" << CGoGNendl;
			return false;
		}

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			int x = imageCoordinates[d].getXCoordinate(), y = imageCoordinates[d].getYCoordinate();
			pixelsLeft(x, y) = pixels(x, y);
		}

		for(int j = height-1; j >= 0; --j)
		{
			for(int i = 0; i < width; ++i)
			{
				out << pixelsLeft(i, j) << " " << std::flush;
			}
			out << std::endl;
		}

		out.close();

		return true;
	}

	return false;
}

bool Surface_DepthMapRendering_Plugin::saveMergedPointCloud(const QString& mapOrigin, const QStringList& mapNames, const QString& directory)
{
	if(!mapOrigin.isEmpty() && !mapNames.empty() && !directory.isEmpty())
	{
		std::vector<PFP2::MAP*> maps;
		maps.reserve(mapNames.size());
		std::vector<VertexAttribute<PFP2::VEC3, PFP2::MAP>> positions;
		for(int i = 0; i < mapNames.size(); ++i)
		{
			MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapNames[i]));
			if(mh_map)
			{
				PFP2::MAP* map = mh_map->getMap();
				maps.push_back(map);

				VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
				if(!position.isValid())
				{
					CGoGNerr << "position attribute is not valid" << CGoGNendl;
					return false;
				}
				positions.push_back(position);
			}
			else
			{
				return false;
			}
		}

		int width = m_depthFBO->getWidth(), height = m_depthFBO->getHeight();

		QString filename(directory);
		filename += "/" + mapOrigin + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += "PointClouds/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += QString::number(width) + "x" + QString::number(height) + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		if(m_correspondance_done)
		{
			filename += mapOrigin + "-" + QString::number(width) + "x" + QString::number(height) + "-Merged-Without.ply";
		}
		else
		{
			filename += mapOrigin + "-" + QString::number(width) + "x" + QString::number(height) + "-Merged-With.ply";
		}

		return Algo::Surface::Export::exportPLYVertMaps<PFP2>(maps, positions, filename.toStdString().c_str(), false);
	}

	return false;
}

void Surface_DepthMapRendering_Plugin::normalEstimation(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin))
	{
		CGoGNout << "Calcul des normales de la carte " << mapGenerated.toStdString() << " .." << CGoGNflush;
		Utils::Chrono chrono;
		chrono.start();

		PFP2::MAP* generated_map = mh_generated->getMap();

		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("position");
		VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_generated->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");

		int width = m_depthFBO->getWidth(), height = m_depthFBO->getHeight();

		Eigen::Matrix<GLfloat, Eigen::Dynamic, 3> matrix;
		matrix.setZero(width*height, 3);

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			matrix.row(imageCoordinates[d].getXCoordinate()*height+imageCoordinates[d].getYCoordinate())
					= Eigen::Vector3f(position[d][0], position[d][1], position[d][2]);
		}

		for(int i = 0; i < matrix.rows()-height; ++i)
		{
			float x1 = matrix(i+height, 0), y1 = matrix(i+height, 1), z1 = matrix(i+height,2);
			float x2 = matrix(i+1, 0), y2 = matrix(i+1, 1), z2 = matrix(i+1,2);
			if((fabs(x1) > FLT_EPSILON || fabs(y1) > FLT_EPSILON || fabs(z1) > FLT_EPSILON)
			   && (fabs(x2) > FLT_EPSILON || fabs(y2) > FLT_EPSILON || fabs(z2) > FLT_EPSILON))
			{
				Eigen::Vector3f u = matrix.row(i+height)-matrix.row(i);	//(x+1;y)-(x;y)
				Eigen::Vector3f v = matrix.row(i+1)-matrix.row(i);	//(x;y+1)-(x;y)
				matrix.row(i) = (u.cross(v)).normalized();
			}
		}

		VertexAttribute<PFP2::VEC3, PFP2::MAP> normal = mh_generated->addAttribute<PFP2::VEC3, VERTEX>("normal");

		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			int pos_in_mat = imageCoordinates[d].getXCoordinate()*height+imageCoordinates[d].getYCoordinate();
			normal[d] = PFP2::VEC3(matrix(pos_in_mat, 0), matrix(pos_in_mat, 1), matrix(pos_in_mat, 2));
			if(normal[d] == position[d])
			{
				generated_map->deleteVertex(d);
			}
		}

		mh_generated->notifyConnectivityModification(false);
		mh_generated->notifyAttributeModification(normal, false);
		mh_generated->updateBB(position);

		CGoGNout << ".. fait en " << chrono.elapsed() << " ms" << CGoGNendl;
	}
}

void Surface_DepthMapRendering_Plugin::confidenceEstimation(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin))
	{
		CGoGNout << "Calcul des valeurs de confiance de visibilité de la carte " << mapGenerated.toStdString() << " .." << CGoGNflush;
		Utils::Chrono chrono;
		chrono.start();

		PFP2::MAP* generated_map = mh_generated->getMap();

		MapParameters& mapParams = m_mapParameterSet[mh_origin];
		Camera* camera = mapParams.depthCameraSet[mapGenerated];

		PFP2::VEC3 position_camera(camera->position().x, camera->position().y, camera->position().z);

		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("position");
		VertexAttribute<PFP2::VEC3, PFP2::MAP> normal = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("normal");
		VertexAttribute<float, PFP2::MAP> visibilityConfidence = mh_generated->addAttribute<float, VERTEX>("VisibilityConfidence");

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			PFP2::VEC3 u = position_camera-position[d];
			u.normalize();
			visibilityConfidence[d] = u*normal[d];
			if(visibilityConfidence[d] != visibilityConfidence[d])
			{	//visibilityConfidence[d]==NaN
				visibilityConfidence[d] = 0.f;
			}
		}

		generated_map->removeAttribute(normal);

		mh_generated->notifyAttributeModification(visibilityConfidence, false);

		CGoGNout << ".. fait en " << chrono.elapsed() << " ms" << CGoGNendl;
	}
}

void Surface_DepthMapRendering_Plugin::findCorrespondingPoints(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin))
	{
		CGoGNout << "Visibilité des points de la carte " << mapGenerated.toStdString() << " .." << CGoGNflush;
		Utils::Chrono chrono;
		chrono.start();

		MapParameters& mapParams = m_mapParameterSet[mh_origin];

		Camera* camera = mapParams.depthCameraSet[mapGenerated];
		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>& existing_depthImage = mapParams.depthImageSet[mapGenerated];

		PFP2::MAP* generated_map = mh_generated->getMap();
		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("position");
		VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_generated->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
		VertexAttribute<float, PFP2::MAP> visibilityConfidence = mh_generated->getAttribute<float, VERTEX>("VisibilityConfidence");

		int width = m_depthFBO->getWidth(), height = m_depthFBO->getHeight();

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> depthImage;
		depthImage.setZero(width, height);

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> maxConfidenceValues;
		maxConfidenceValues.setZero(width, height);

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> confidenceValues;
		confidenceValues.setZero(width, height);

		//Change current rendering camera, and reposition it correctly
		qglviewer::Vec camera_position = camera->position();
		m_schnapps->getSelectedView()->setCurrentCamera(camera, false);
		camera->setPosition(camera_position);

		GLfloat n = fabs(camera->zNear()), f = fabs(camera->zFar());

		/*
		 * z_w = valeur de profondeur de l'image ; dans l'intervalle [0;1]
		 * n = distance non signée du zNear à la caméra
		 * f = distance non signée du zFar à la caméra
		 * p = (f*n)/(z_w*(f-n)-f) => Fonction qui calcule la distance d'un point à la caméra, selon la valeur de profondeur de la carte de profondeur
		 * - ((f*n) * (f-n)) / ((z_w*(f-n)-f)*(z_w*(f-n)-f)) => Dérivée 1ère de la fonction précédente
		 * z_w = (f*(p+n))/(p*(f-n)) => Fonction inverse qui pour une distance d'un point à la caméra, donne sa valeur de profondeur
		 */

		for(QHash<QString, Camera*>::iterator it = mapParams.depthCameraSet.begin(); it != mapParams.depthCameraSet.end(); ++it)
		{
			if(it.key().compare(mh_generated->getName()) != 0)
			{
				MapHandler<PFP2>* mh_current = static_cast<MapHandler<PFP2>*>(mapParams.projectedMapSet[it.key()]);

				m_shaderScalarFieldReal->setAttributePosition(mh_current->getVBO("position"));
				m_shaderScalarFieldReal->setAttributeScalar(mh_current->getVBO("VisibilityConfidence"));

				m_depthFBO->bind();
				glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );	//To clean the color and depth textures

				mh_current->draw(m_shaderScalarFieldReal, CGoGN::Algo::Render::GL2::POINTS);

				glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, depthImage.data());	//Read pixels from depth map
				glReadPixels(0, 0, width, height, GL_RED, GL_FLOAT, confidenceValues.data()); //Read pixels from color map to get confidence values
				m_depthFBO->unbind();

				depthImage = depthImage.array()*2-1;    //Set range to [-1;1]

				depthImage = (existing_depthImage.array()-depthImage.array()).abs();

				for(int i = 0; i < depthImage.rows(); ++i)
				{
					for(int j = 0; j < depthImage.cols(); ++j)
					{
						GLfloat z_w = 0.5*existing_depthImage(i, j)+0.5;	//Set value in [0;1]
						GLdouble upper_bound = (f*n)/(z_w*(f-n)-f);
						GLdouble lower_bound = upper_bound;
						upper_bound += (f-n)/100.f, lower_bound -= (f-n)/100.f;
						upper_bound = (f*((upper_bound)+n))/((upper_bound)*(f-n)), lower_bound = (f*((lower_bound)+n))/((lower_bound)*(f-n));
						if(depthImage(i, j) < fabs(upper_bound-lower_bound))
						{
							maxConfidenceValues(i, j) = std::max(maxConfidenceValues(i, j), confidenceValues(i, j));
						}
					}
				}
			}
		}

		m_schnapps->getSelectedView()->setCurrentCamera("camera_0", false);

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			int x = imageCoordinates[d].getXCoordinate(), y = imageCoordinates[d].getYCoordinate();
			if(visibilityConfidence[d] < maxConfidenceValues(x, y))
			{
				generated_map->deleteVertex(d);
			}
		}

		mh_generated->notifyConnectivityModification(false);
		mh_generated->notifyAttributeModification(position, false);
		mh_generated->notifyAttributeModification(visibilityConfidence, false);
		mh_generated->notifyAttributeModification(imageCoordinates, false);
		mh_generated->updateBB(position);

		CGoGNout << ".. fait en " << chrono.elapsed() << " ms" << CGoGNendl;

		m_correspondance_done = true;
	}
}

void Surface_DepthMapRendering_Plugin::deleteBackground(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin))
	{
		MapParameters& mapParams = m_mapParameterSet[mh_origin];
		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>& existing_depthImage = mapParams.depthImageSet[mapGenerated];

		PFP2::MAP* generated_map = mh_generated->getMap();

		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("position");
		VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_generated->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");

		TraversorF<PFP2::MAP> trav_face_map(*generated_map);
		Dart next;
		for(Dart d = trav_face_map.begin(); d != trav_face_map.end(); d = next)
		{
			next = trav_face_map.next();
			bool stop = false;
			Traversor2FV<PFP2::MAP> trav_vert_face_map(*generated_map, d);
			for(Dart dd = trav_vert_face_map.begin(); !stop && dd != trav_vert_face_map.end(); dd = trav_vert_face_map.next())
			{
				GLfloat color = existing_depthImage(imageCoordinates[dd].getXCoordinate(),imageCoordinates[dd].getYCoordinate());
				if(fabs(1 - color) < FLT_EPSILON)
				{
					//Le point fait partie du fond de l'image
					//position[d] = PFP2::VEC3(0.f, 0.f, 0.f);
					generated_map->deleteFace(d);
					stop = true;
				}
			}
		}

		mh_generated->notifyConnectivityModification(false);
		mh_generated->notifyAttributeModification(position, false);
		mh_generated->notifyAttributeModification(imageCoordinates, false);
		mh_generated->updateBB(position);
	}
}

void Surface_DepthMapRendering_Plugin::exportModelPly(const QString& mapName, const QString& directory)
{
	MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
	MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);

	if(mh_map && !mapName.isEmpty() && !directory.isEmpty())
	{
		PFP2::MAP* map = mh_map->getMap();
		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");

		QString filename(directory);
		filename += "/" + mapName + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += mapName + ".ply";

		Algo::Surface::Export::exportPLY<PFP2>(*map, position, filename.toStdString().c_str(), false);
	}
}

bool Surface_DepthMapRendering_Plugin::moveDownDecomposition(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin))
	{
		MapParameters& mapParams = m_mapParameterSet[mh_origin];
	}
	return false;
}

#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_Plugin, Surface_DepthMapRendering_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_PluginD, Surface_DepthMapRendering_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
