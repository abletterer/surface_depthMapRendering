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

	m_shaderSimpleColor = new CGoGN::Utils::ShaderSimpleColor();
	m_shaderScalarFieldReal = new CGoGN::Utils::ShaderScalarFieldReal();

	registerShader(m_shaderSimpleColor);
	registerShader(m_shaderScalarFieldReal);

	connect(m_depthMapRenderingAction, SIGNAL(triggered()), this, SLOT(openDepthMapRenderingDialog()));

	connect(m_depthMapRenderingDialog, SIGNAL(accepted()), this, SLOT(closeDepthMapRenderingDialog()));
	connect(m_depthMapRenderingDialog->button_cancel, SIGNAL(clicked()), this, SLOT(closeDepthMapRenderingDialog()));
	connect(m_depthMapRenderingDialog->button_ok, SIGNAL(clicked()), this, SLOT(closeDepthMapRenderingDialog()));
	connect(m_depthMapRenderingDialog->button_moveDown, SIGNAL(clicked()), this, SLOT(moveDownFromDialog()));
	connect(m_depthMapRenderingDialog->button_moveUp, SIGNAL(clicked()), this, SLOT(moveUpFromDialog()));

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

	disconnect(m_depthMapRenderingDialog, SIGNAL(accepted()), this, SLOT(closeDepthMapRenderingDialog()));
	disconnect(m_depthMapRenderingDialog->button_cancel, SIGNAL(clicked()), this, SLOT(closeDepthMapRenderingDialog()));
	disconnect(m_depthMapRenderingDialog->button_ok, SIGNAL(clicked()), this, SLOT(closeDepthMapRenderingDialog()));
	disconnect(m_depthMapRenderingDialog->button_moveDown, SIGNAL(clicked()), this, SLOT(moveDownFromDialog()));
	disconnect(m_depthMapRenderingDialog->button_moveUp, SIGNAL(clicked()), this, SLOT(moveUpFromDialog()));

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

void Surface_DepthMapRendering_Plugin::render(const QString& mapName, bool saveData, const QString& directory)
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

		std::vector<QString> mapNames;
		mapNames.reserve(mapParams.depthCameraSet.size());

		for(QHash<QString, Camera*>::iterator it = mapParams.depthCameraSet.begin(); it != mapParams.depthCameraSet.end(); ++it)
		{
			Camera* camera = it.value();
			QString cameraName(camera->getName());

			QString generatedName(mapName);
			generatedName += "-" + cameraName;

			mapNames.push_back(generatedName);

			m_schnapps->getSelectedView()->setCurrentCamera(camera, false);

			m_depthFBO->bind();
			glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );	//To clean the color and depth textures
			glBindTexture(GL_TEXTURE_2D, *m_depthFBO->getDepthTexId());

			mh_map->draw(m_shaderSimpleColor, CGoGN::Algo::Render::GL2::TRIANGLES);	//Render the map into the FrameBufferObject

			//Read pixels of the generated texture and store them in an array
			glGetTexImage(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, GL_FLOAT, pixels.data());
			glBindTexture(GL_TEXTURE_2D, 0);
			m_depthFBO->unbind();

			m_schnapps->getSelectedView()->setCurrentCamera("camera_0", false);

			MapHandlerGen* mhg_generated = m_schnapps->addMap(generatedName, 2);
			mapParams.projectedMapSet[generatedName] = mhg_generated;

			mapParams.decompositionLevel[generatedName] = 0;

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

			pixels = pixels.array()*2-1;	//Put depth values in the range [-1;1]

			mapParams.depthImageSet[generatedName] = pixels;

			generated_map->enableQuickTraversal<PFP2::MAP, VERTEX>();

			mh_generated->notifyConnectivityModification(false);
			mh_generated->notifyAttributeModification(planeCoordinatesGenerated, false);
			mh_generated->notifyAttributeModification(imageCoordinatesGenerated, false);
			project2DImageTo3DSpace(mapName, generatedName);
		}

		CGoGNout << "Temps d'échantillonnage : " << chrono.elapsed() << " ms " << CGoGNflush;
		CGoGNout << "pour " << mapParams.depthCameraSet.size() << " vue(s) différente(s) " << CGoGNflush;
		CGoGNout << "de taille " << width << "x" << height << CGoGNflush;
		CGoGNout << " sur un objet composé de " << mh_map->getMap()->getNbCells(VERTEX) << " point(s)" << CGoGNendl;

		if(saveData)
		{
			saveMergedPointCloud(mapName, mapNames, directory);
		}

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

		VertexAttribute<PFP2::VEC3, PFP2::MAP> positionGenerated = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("position");
		if(!positionGenerated.isValid())
		{
			positionGenerated = mh_generated->addAttribute<PFP2::VEC3, VERTEX>("position");
		}

		VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinatesGenerated = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
		if(!planeCoordinatesGenerated.isValid())
		{
			CGoGNerr << "PlaneCoordinates attribute is not valid" << CGoGNendl;
			return;
		}

		VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinatesGenerated = mh_generated->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
		if(!imageCoordinatesGenerated.isValid())
		{
			CGoGNerr << "ImageCoordinates attribute is not valid" << CGoGNendl;
			return;
		}

		Camera* camera = mapParams.depthCameraSet[mapGenerated];
		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> pixels = mapParams.depthImageSet[mapGenerated];

		TraversorF<PFP2::MAP> trav_face_map(*generated_map);
		Dart next;
		for(Dart d = trav_face_map.begin(); d != trav_face_map.end(); d = next)
		{
			next = trav_face_map.next();
			bool stop = false;
			Traversor2FV<PFP2::MAP> trav_vert_face_map(*generated_map, d);
			for(Dart dd = trav_vert_face_map.begin(); !stop && dd != trav_vert_face_map.end(); dd = trav_vert_face_map.next())
			{
				float color = pixels(imageCoordinatesGenerated[dd].getXCoordinate(),imageCoordinatesGenerated[dd].getYCoordinate());
				if(fabs(1-color)<FLT_EPSILON)
				{
					//Le point fait partie du fond de l'image
					generated_map->deleteFace(d);
					stop = true;
				}
			}
		}
		generated_map->updateQuickTraversal<PFP2::MAP, VERTEX>();

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
			float color = pixels(imageCoordinatesGenerated[d].getXCoordinate(),imageCoordinatesGenerated[d].getYCoordinate());

			PFP2::VEC4 pos = PFP2::VEC4(planeCoordinatesGenerated[d][0], planeCoordinatesGenerated[d][1], color, 1.f);

			pos = model_view_projection_matrix_inv*pos;

			positionGenerated[d] = PFP2::VEC3(pos[0]/pos[3], pos[1]/pos[3], pos[2]/pos[3]);
		}

		mh_generated->notifyAttributeModification(positionGenerated, false);
		mh_generated->updateBB(positionGenerated);
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

		filename += mapGenerated + ".ply";

		return Algo::Surface::Export::exportPLYVert<PFP2>(*generated_map, position, filename.toStdString().c_str(), false);
	}

	return false;
}

bool Surface_DepthMapRendering_Plugin::saveDepthMap(const QString& mapOrigin, const QString& mapGenerated, const QString& directory)
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
		out.open(filename.toStdString() + "-depthMap.dat", std::ios::out);
		if(!out.good())
		{
			CGoGNerr << "Unable to open file" << CGoGNendl;
			return false;
		}

		for(unsigned int j = m_depthFBO->getHeight()-1; j >= 0; --j)
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

bool Surface_DepthMapRendering_Plugin::saveMergedPointCloud(const QString& mapOrigin, const std::vector<QString>& mapNames, const QString& directory)
{
	if(!mapOrigin.isEmpty() && !mapNames.empty() && !directory.isEmpty())
	{
		std::vector<PFP2::MAP*> maps;
		maps.reserve(mapNames.size());
		std::vector<VertexAttribute<PFP2::VEC3, PFP2::MAP>> positions;
		for(unsigned int i = 0; i < mapNames.size(); ++i)
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

		QString filename(directory);
		filename += "/" + mapOrigin + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += "PointClouds/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += QString::number(m_depthFBO->getWidth()) + "x" + QString::number(m_depthFBO->getHeight()) + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += mapOrigin + "-Merged.ply";

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

		for(int i = 0; i < matrix.rows()-1; ++i)
		{
			Eigen::Vector3f u = matrix.row(i+height)-matrix.row(i);	//(x+1;y)-(x;y)
			Eigen::Vector3f v = matrix.row(i+1)-matrix.row(i);	//(x;y+1)-(x;y)
			matrix.row(i) = (u.cross(v)).normalized();
		}

		VertexAttribute<PFP2::VEC3, PFP2::MAP> normal = mh_generated->addAttribute<PFP2::VEC3, VERTEX>("normal");

		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			int pos_in_mat = imageCoordinates[d].getXCoordinate()*height+imageCoordinates[d].getYCoordinate();
			normal[d] = PFP2::VEC3(matrix(pos_in_mat, 0), matrix(pos_in_mat, 1), matrix(pos_in_mat, 2));
		}

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

		int width = m_depthFBO->getWidth(), height = m_depthFBO->getHeight();

		Eigen::Matrix<GLfloat, Eigen::Dynamic, 3> matrix;
		matrix.setZero(width*height, 3);

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			visibilityConfidence[d] = (position_camera-position[d])*(normal[d]);
			if(visibilityConfidence[d] != visibilityConfidence[d])
			{	//visibilityConfidence[d]==NaN
				visibilityConfidence[d] = 0.f;
			}
		}

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
		CGoGNout << "Appariement des points de la carte " << mapGenerated.toStdString() << " :" << CGoGNendl;

		MapParameters& mapParams = m_mapParameterSet[mh_origin];

		int width = m_depthFBO->getWidth(), height = m_depthFBO->getHeight();

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> depthImage;
		depthImage.setZero(width, height);

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> confidenceValues;
		confidenceValues.setZero(width, height);

		m_shaderScalarFieldReal->setAttributePosition(mhg_generated->getVBO("position"));
		m_shaderScalarFieldReal->setAttributeScalar(mhg_generated->getVBO("VisibilityConfidence"));

		for(QHash<QString, Camera*>::iterator it = mapParams.depthCameraSet.begin(); it != mapParams.depthCameraSet.end(); ++it)
		{
			if(it.key().compare(mh_generated->getName()) != 0)
			{
				Camera* current_camera = it.value();
				CGoGNout << "\t Par rapport à la caméra " << current_camera->getName().toStdString() << " .." << CGoGNflush;
				Utils::Chrono chrono;
				chrono.start();

				MapHandler<PFP2>* mh_current = static_cast<MapHandler<PFP2>*>(mapParams.projectedMapSet[it.key()]);

				m_schnapps->getSelectedView()->setCurrentCamera(current_camera, false);

				m_depthFBO->bind();
				glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );	//To clean the color and depth textures

				mh_generated->draw(m_shaderScalarFieldReal, CGoGN::Algo::Render::GL2::POINTS);	//Render the map into the FrameBufferObject

				//Read pixels and store them in an array
				glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, depthImage.data());
				glReadPixels(0, 0, width, height, GL_RED, GL_FLOAT, confidenceValues.data());
				m_depthFBO->unbind();

				m_schnapps->getSelectedView()->setCurrentCamera("camera_0", false);

				depthImage = depthImage.array()*2-1;

				Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> current_depthImage = mapParams.depthImageSet[mh_current->getName()];

				for(int i = 0; i < confidenceValues.rows(); ++i)
				{
					for(int j = 0; j < confidenceValues.cols(); ++j)
					{
						float value = fabs(depthImage(i,j)-current_depthImage(i,j));
					}
				}

				CGoGNout << ".. fait en " << chrono.elapsed() << " ms" << CGoGNendl;
			}
		}
	}
}

//void Surface_DepthMapRendering_Plugin::applyFilter(const QString& mapOrigin, const QString& mapGenerated, std::vector<int>& f, int size_x, int size_y)
//{
//	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
//	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

//	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
//	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

//	if(mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin) && !f.empty())
//	{
//		MapParameters& mapParams = m_mapParameterSet[mh_origin];

//		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>& pixels = mapParams.depthImageSet[mapGenerated];
//		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> res(pixels);

//		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> filter;
//		filter.setZero(size_x, size_y);

//		for(int i = 0; i < filter.rows(); ++i)
//		{
//			for(int j = 0; j < filter.cols(); ++j)
//			{
//				filter(i, j) = f[j+size_x*i];
//			}
//		}

//		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> tmp;
//		tmp.setZero(size_x, size_y);


//		for(int i = 1; i < pixels.rows()-1; ++i)
//		{
//			for(int j = 1; j < pixels.cols()-1; ++j)
//			{
//				tmp <<	pixels(i-1, j-1), pixels(i-1, j), pixels(i-1, j+1),
//						pixels(i, j-1), pixels(i, j), pixels(i, j+1),
//						pixels(i+1, j-1), pixels(i+1, j), pixels(i+1, j+1);

//				tmp.array() *= filter.array();

//				res(i, j) = tmp.sum();
//			}
//		}
//	}
//}

#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_Plugin, Surface_DepthMapRendering_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_PluginD, Surface_DepthMapRendering_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
