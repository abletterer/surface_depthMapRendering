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

	connect(m_depthMapRenderingDialog->button_lower_resolution, SIGNAL(clicked()), this, SLOT(lowerResolutionFromDialog()));
	connect(m_depthMapRenderingDialog->button_upper_resolution, SIGNAL(clicked()), this, SLOT(upperResolutionFromDialog()));
	connect(m_depthMapRenderingDialog->button_saveMergedPointCloud, SIGNAL(clicked()), this, SLOT(saveMergedPointCloudFromDialog()));

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
	if(m_fbo)
	{
		delete m_fbo;
	}

	disconnect(m_depthMapRenderingAction, SIGNAL(triggered()), this, SLOT(openDepthMapRenderingDialog()));

	disconnect(m_depthMapRenderingDialog->button_lower_resolution, SIGNAL(clicked()), this, SLOT(lowerResolutionFromDialog()));
	disconnect(m_depthMapRenderingDialog->button_upper_resolution, SIGNAL(clicked()), this, SLOT(upperResolutionFromDialog()));
	disconnect(m_depthMapRenderingDialog->button_saveMergedPointCloud, SIGNAL(clicked()), this, SLOT(saveMergedPointCloudFromDialog()));

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

void Surface_DepthMapRendering_Plugin::lowerResolutionFromDialog()
{
	QList<QListWidgetItem*> currentItems = m_depthMapRenderingDialog->list_maps->selectedItems();
	if(!currentItems.empty())
	{
		QString mapName = currentItems[0]->text();
		MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
		if(mhg_map && m_mapParameterSet.contains(mhg_map))
		{
			MapParameters& mapParams = m_mapParameterSet[mhg_map];
			for(QHash<QString, MapHandlerGen*>::iterator it = mapParams.projectedMapSet.begin(); it != mapParams.projectedMapSet.end(); ++it)
			{
				lowerResolution(mapName, it.key());
				deleteBackground(mapName, it.key());
			}
			m_schnapps->getSelectedView()->updateGL();
		}
		else
		{
			CGoGNout << "Nothing to do for this map" << CGoGNendl;
		}
	}
}

void Surface_DepthMapRendering_Plugin::upperResolutionFromDialog()
{
	QList<QListWidgetItem*> currentItems = m_depthMapRenderingDialog->list_maps->selectedItems();
	if(!currentItems.empty())
	{
		QString mapName = currentItems[0]->text();
		MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
		if(mhg_map && m_mapParameterSet.contains(mhg_map))
		{
			MapParameters& mapParams = m_mapParameterSet[mhg_map];
			for(QHash<QString, int>::iterator it = mapParams.decompositionLevelSet.begin(); it != mapParams.decompositionLevelSet.end(); ++it)
			{
				upperResolution(mapName, it.key());
				deleteBackground(mapName, it.key());
				m_schnapps->getSelectedView()->updateGL();
			}
		}
		else
		{
			CGoGNout << "Nothing to do for this map" << CGoGNendl;
		}
	}
}

void Surface_DepthMapRendering_Plugin::saveMergedPointCloudFromDialog()
{
	QList<QListWidgetItem*> currentItems = m_depthMapRenderingDialog->list_maps->selectedItems();
	if(!currentItems.empty())
	{
		QString mapName = currentItems[0]->text();
		MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
		if(mhg_map && m_mapParameterSet.contains(mhg_map))
		{
			MapParameters& mapParams = m_mapParameterSet[mhg_map];
			QStringList mapNames = mapParams.depthCameraSet.keys();
			saveMergedPointCloud(mapName, mapNames);
		}
	}
}

void Surface_DepthMapRendering_Plugin::mapAdded(MapHandlerGen *map)
{
	connect(map, SIGNAL(vboRemoved(Utils::VBO*)), this, SLOT(vboRemoved(Utils::VBO*)));
	connect(map, SIGNAL(selectedCellsChanged(CellSelectorGen*)), this, SLOT(selectedCellsChanged(CellSelectorGen*)));
}

void Surface_DepthMapRendering_Plugin::mapRemoved(MapHandlerGen *map)
{
	disconnect(map, SIGNAL(vboRemoved(Utils::VBO*)), this, SLOT(vboRemoved(Utils::VBO*)));
	disconnect(map, SIGNAL(selectedCellsChanged(CellSelectorGen*)), this, SLOT(selectedCellsChanged(CellSelectorGen*)));
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

void Surface_DepthMapRendering_Plugin::selectedCellsChanged(CellSelectorGen* cs)
{
	MapHandlerGen* mhg_map = static_cast<MapHandlerGen*>(QObject::sender());
	MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);

	if(mh_map && m_main_object)
	{
		VertexAttribute<NoTypeNameAttribute<std::vector<PointCorrespondance>>, PFP2::MAP> correspondingPointsAttribute
				= mh_map->getAttribute<NoTypeNameAttribute<std::vector<PointCorrespondance>>, VERTEX>("CorrespondingPoints");
		if(!correspondingPointsAttribute.isValid())
		{
			CGoGNerr << "CorrespondingPoints attribute does not exist" << CGoGNendl;
			return;
		}

		CellSelector<PFP2::MAP, VERTEX>* vertexSelector = static_cast<CellSelector<PFP2::MAP, VERTEX>*>(cs);
		if(vertexSelector->getNbSelectedCells() == 0)
		{
			//Démarquage de tous les sommets de toutes les cartes
			MapParameters& mapParams = m_mapParameterSet[m_main_object];
			for(QHash<QString, MapHandlerGen*>::iterator it = mapParams.projectedMapSet.begin(); it != mapParams.projectedMapSet.end(); ++it)
			{
				if(it.value() != mh_map)
				{
					MapHandler<PFP2>* mh_current = static_cast<MapHandler<PFP2>*>(it.value());
					CellSelector<PFP2::MAP, VERTEX>* cur_selector = mh_current->getCellSelector<VERTEX>("selector");
					if(cur_selector->getNbSelectedCells() > 0)
					{
						cur_selector->unselect(cur_selector->getSelectedCells());
					}
				}
			}
		}
		else
		{
			Vertex selected_point = vertexSelector->getSelectedCells()[0];

			std::vector<PointCorrespondance>& correspondingPoints = correspondingPointsAttribute[selected_point];
			for(std::vector<PointCorrespondance>::const_iterator it = correspondingPoints.begin(); it!=correspondingPoints.end(); ++it)
			{
				PointCorrespondance tmp = *it;
				MapHandler<PFP2>* mh_current = static_cast<MapHandler<PFP2>*>(tmp.map);

				CellSelector<PFP2::MAP, VERTEX>* cur_selector = mh_current->getCellSelector<VERTEX>("selector");
				if(cur_selector->getNbSelectedCells() > 0)
				{
					cur_selector->unselect(cur_selector->getSelectedCells());
				}
				cur_selector->select(tmp.vertex, false);
			}
		}

		m_schnapps->getSelectedView()->updateGL();
	}
}

void Surface_DepthMapRendering_Plugin::createFBO(int width, int height)
{
	m_fbo = new CGoGN::Utils::FBO(width, height);
	m_fbo->createAttachDepthTexture();
	m_fbo->createAttachColorTexture(GL_R32F);
}

void Surface_DepthMapRendering_Plugin::changePositionVBO(const QString& view, const QString& map, const QString& vbo)
{
	View* v = m_schnapps->getView(view);
	MapHandlerGen* m = m_schnapps->getMap(map);
	if(v && m)
	{
		m_main_object = m;
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
		positions.push_back(qglviewer::Vec(0,-1,2)); // Swap
		positions.push_back(qglviewer::Vec(0,1,-2));	// Swap
		positions.push_back(qglviewer::Vec(0,-1,-2));

		positions.push_back(qglviewer::Vec(1,2,0));
		positions.push_back(qglviewer::Vec(1,-2,0));
		positions.push_back(qglviewer::Vec(-1,2,0));
		positions.push_back(qglviewer::Vec(-1,-2,0));

		positions.push_back(qglviewer::Vec(2,0,1));
		positions.push_back(qglviewer::Vec(2,0,-1));
		positions.push_back(qglviewer::Vec(-2,0,1));
		positions.push_back(qglviewer::Vec(-2,0,-1));

		qglviewer::Vec bb_min, bb_max;

		mh_map->transformedBB(bb_min, bb_max);

		qglviewer::Vec center = (bb_min+bb_max)/2.f;

		for(int i = 0; i < nbMax; ++i)
		{
			QString cameraName(baseName);
			cameraName.append(QString::number(i));
			Camera* camera = m_schnapps->addCamera(cameraName);

			qglviewer::Vec camera_position(camera->position());

			float radius = 1;

			radius *= (bb_max-bb_min).norm()/4;

			camera_position.x = center.x + radius*positions[i].x;
			camera_position.y = center.y + radius*positions[i].y;
			camera_position.z = center.z + radius*positions[i].z;

			camera->setPosition(camera_position);

			camera->lookAt(center);

			camera->setSceneCenter(center);
			camera->centerScene();

			camera_position = camera->position();

			camera->disableViewsBoundingBoxFitting();

			camera->setZNear((camera_position-center-radius*positions[i]/2).norm());
			camera->setZFar((camera_position-center+radius*positions[i]/2).norm());
			camera->setStandard(false);

			QString generatedName(mapName);
			generatedName += "-" + cameraName;

			mapParams.depthCameraSet[generatedName] = camera;
		}
	}
}

void Surface_DepthMapRendering_Plugin::render(const QString& mapName, const QString& directory)
{
	MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
	MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);

	if(m_fbo && mh_map && m_mapParameterSet.contains(mhg_map))
	{
		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
		if(!position.isValid())
		{
			CGoGNerr << "position attribute is not valid" << CGoGNendl;
			return;
		}

		MapParameters& mapParams = m_mapParameterSet[mhg_map];

		const int width = m_fbo->getWidth(), height = m_fbo->getHeight();
		m_shaderSimpleColor->setAttributePosition(mapParams.positionVBO);

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> pixels(width, height);

		int total_sampling_time = 0, total_saving_time = 0;

		Camera* o_camera = m_schnapps->getSelectedView()->getCurrentCamera();

		for(QHash<QString, Camera*>::iterator it = mapParams.depthCameraSet.begin(); it != mapParams.depthCameraSet.end(); ++it)
		{
			Camera* camera = it.value();
			QString cameraName(camera->getName());

			QString generatedName(mapName);
			generatedName += "-" + cameraName;

			qglviewer::Vec camera_position = camera->position();

			m_schnapps->getSelectedView()->setCurrentCamera(camera);

			Utils::Chrono chrono;
			chrono.start();

			m_fbo->bind();
			glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );	//To clean the color and depth textures

			mh_map->draw(m_shaderSimpleColor, CGoGN::Algo::Render::GL2::TRIANGLES);	//Render the map into the FrameBufferObject

			glBindTexture(GL_TEXTURE_2D, *m_fbo->getDepthTexId());
			glGetTexImage(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, GL_FLOAT, pixels.data());
			m_fbo->unbind();
			
			total_sampling_time += chrono.elapsed();
			
            QString filename(directory);
			filename += "/" + mapName + "/";
			mkdir(filename.toStdString().c_str(), 0777);

			filename += "DepthMaps/";
			mkdir(filename.toStdString().c_str(), 0777);

			filename += QString::number(width) + "x" + QString::number(height) + "/";
			mkdir(filename.toStdString().c_str(), 0777);

			filename += generatedName;
			
			chrono.start();
	
			std::ofstream out;
			out.open(filename.toStdString() + "-originalDepthMap.dat", std::ios::out);
			if(!out.good())
			{
				CGoGNerr << "Unable to open file" << CGoGNendl;
				return;
			}

			out << pixels;

			out.close();

			out.open(filename.toStdString() + "-MVPMatrix.dat", std::ios::out);
			if(!out.good())
			{
				CGoGNerr << "Unable to open file" << CGoGNendl;
				return;
			}

			GLdouble mvp_matrix[16];
			camera->getModelViewProjectionMatrix(mvp_matrix);

			for(int i = 0; i < 4; ++i)
			{
				for(int j = 0; j < 4; ++j)
				{
					out << mvp_matrix[i+j*4] << " " << std::flush;
				}
				out << std::endl;
			}

			out.close();
			
			total_saving_time += chrono.elapsed();

//			m_schnapps->getSelectedView()->setCurrentCamera("camera_0");

//			MapHandlerGen* mhg_generated = m_schnapps->addMap(generatedName, 2);
//			mapParams.projectedMapSet[generatedName] = mhg_generated;

////			pixels.array() = pixels.array()*2-1;	//Put depth values in the range [-1;1]

//			mapParams.depthImageSet[generatedName] = pixels;
//			mapParams.decompositionLevelSet[generatedName] = 0;

//			MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);
//			PFP2::MAP* generated_map = mh_generated->getMap();

//			VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinatesGenerated =
//				mh_generated->addAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
//			VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinatesGenerated =
//				mh_generated->addAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");

//			Algo::Surface::Tilings::Square::Grid<PFP2> grid(*generated_map, width-1, height-1);
//			grid.embedIntoGrid(planeCoordinatesGenerated, 2, 2);

//			std::vector<Dart>& vDarts = grid.getVertexDarts();

//			for(int i = 0; i < width; ++i)
//			{
//				for(int j = 0; j < height; ++j)
//				{
//					imageCoordinatesGenerated[vDarts[j*width+i]].setCoordinates(i, j);
//				}
//			}
//			project2DImageTo3DSpace(mapName, generatedName);
//			deleteBackground(mapName, generatedName);

//			mh_generated->notifyAttributeModification(planeCoordinatesGenerated, false);
//			mh_generated->notifyAttributeModification(imageCoordinatesGenerated, false);
		}

//		m_schnapps->getSelectedView()->setCurrentCamera(o_camera);

		CGoGNout << "Temps d'échantillonnage : " << total_sampling_time << " ms " << CGoGNflush;
		CGoGNout << "pour " << mapParams.depthCameraSet.size() << " vue(s) différente(s) " << CGoGNflush;
		CGoGNout << "de taille " << width << "x" << height << CGoGNflush;
		CGoGNout << " sur un objet composé de " << mh_map->getMap()->getNbCells(VERTEX) << " point(s)" << CGoGNendl;
		CGoGNout << "Temps d'enregistrement des cartes de profondeur : " << total_saving_time << " ms " << CGoGNendl;

//		m_schnapps->getSelectedView()->updateGL();
	}
}

void Surface_DepthMapRendering_Plugin::project2DImageTo3DSpace(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);

	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(m_fbo && mh_origin && mh_generated && m_mapParameterSet.contains(mhg_origin))
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

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>& pixels = mapParams.depthImageSet[mapGenerated];

		GLdouble mvp_matrix[16];
		camera->getModelViewProjectionMatrix(mvp_matrix);

		PFP2::MATRIX44 model_view_projection_matrix, model_view_projection_matrix_inv;

		for(int i = 0; i < 4; ++i)
		{
			for(int j = 0; j < 4; ++j)
			{
				model_view_projection_matrix(i, j) = mvp_matrix[i+j*4];
			}
		}

//		std::cout << "-----" << std::endl;
//		std::cout << camera->getName().toStdString() << std::endl;
//		std::cout << model_view_projection_matrix << std::endl;
//		std::cout << "-----" << std::endl;t

		model_view_projection_matrix.invert(model_view_projection_matrix_inv);

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			float color = pixels(imageCoordinates[d].getXCoordinate(), imageCoordinates[d].getYCoordinate());

			if(fabs(1-color) > FLT_EPSILON)
			{
				PFP2::VEC4 pos = PFP2::VEC4(planeCoordinates[d][0], planeCoordinates[d][1], color, 1.f);

				pos = model_view_projection_matrix_inv*pos;

				position[d] = PFP2::VEC3(pos[0]/pos[3], pos[1]/pos[3], pos[2]/pos[3]);
			}
		}

		mh_generated->notifyConnectivityModification(false);
		mh_generated->notifyAttributeModification(position, false);
	}
}

bool Surface_DepthMapRendering_Plugin::lowerResolution(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mhg_origin))
	{
		MapParameters& mapParams = m_mapParameterSet[mhg_origin];
		int& level = mapParams.decompositionLevelSet[mapGenerated];
		if(pow(2, level+1) < m_fbo->getWidth() && pow(2, level+1) < m_fbo->getHeight())
		{
			Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>& pixels = mapParams.depthImageSet[mapGenerated];

			const int img_width = m_fbo->getWidth(), img_height = m_fbo->getHeight();
			++level;
			const int l_p = pow(2, level);
			const int width = img_width/l_p, height = img_height/l_p;
			const int width2 = width*2, height2 = height*2;

			Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> pixels_tmp(pixels);

			#pragma omp parallel for
			for(int j = 0; j < height2; ++j)
			{
				for(int i = 0; i < width2-1; i += 2)
				{
					pixels(i/2, j) = pixels_tmp(i, j);
				}
				for(int i = 1; i < width2-1; i += 2)
				{
					float impair = pixels_tmp(i, j);
					float pair_1 = pixels_tmp(i-1, j);
					float pair_2 = pixels_tmp(i+1, j);

					impair -= (pair_1+pair_2)/2.f;
					pixels(width+i/2, j) = impair/2.f;
				}
			}

			//Traitement spécifique pour la dernière colonne (différence avec le pair situé à gauche)
			#pragma omp parallel for
			for(int j = 0; j < height2; ++j)
			{
				pixels(width2-1, j) = (pixels(width2-1, j)-pixels_tmp(width2-2, j))/2.f;
			}

			pixels_tmp = Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>(pixels);

			#pragma omp parallel for
			for(int i = 0; i < width2; ++i)
			{
				for(int j = 0; j < height2-1; j += 2)
				{
					pixels(i, j/2) = pixels_tmp(i, j);
				}
				for(int j = 1; j < height2-1; j += 2)
				{
					float impair = pixels_tmp(i, j);
					float pair_1 = pixels_tmp(i, j-1);
					float pair_2 = pixels_tmp(i, j+1);

					impair -= (pair_1+pair_2)/2.f;
					pixels(i, height+j/2) = impair/2.f;
				}

				//Traitement spécifique pour la dernière ligne (différence avec le pair situé au dessus)
				pixels(i, height2-1) = (pixels(i, height2-1)-pixels_tmp(i, height2-2))/2.f;
			}

			regenerateMap(mapOrigin, mapGenerated);
			project2DImageTo3DSpace(mapOrigin, mapGenerated);
//			mh_generated->updateBB();

			return true;
		}
	}
	return false;
}

bool Surface_DepthMapRendering_Plugin::upperResolution(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mhg_origin))
	{
		MapParameters& mapParams = m_mapParameterSet[mhg_origin];
		int& level = mapParams.decompositionLevelSet[mapGenerated];

		if(level>0)
		{
			Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>& pixels = mapParams.depthImageSet[mapGenerated];

			const int img_width = m_fbo->getWidth(), img_height = m_fbo->getHeight();
			--level;
			const int l_p = pow(2, level);
			const int width = img_width/l_p, height = img_height/l_p;

			Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> pixels_tmp(pixels);


			#pragma omp parallel for
			for(int j = 0; j < height; ++j)
			{
				for(int i = 0; i < width-1; i += 2)
				{
					pixels(i, j) = pixels_tmp(i/2, j);
				}

				for(int i = 1, i_2 = 0; i < width-1; i += 2, ++i_2)
				{
					float impair = pixels_tmp(width/2+i_2, j)*2;
					float pair_1 = pixels_tmp(i_2, j);
					float pair_2 = pixels_tmp(i_2+1, j);

					impair += (pair_1+pair_2)/2.f;
					pixels(i, j) = impair;
				}
			}

			//Traitement spécifique pour la dernière colonne (différence avec le pair situé à gauche)
			#pragma omp parallel for
			for(int j = 0; j < height; ++j)
			{
				pixels(width-1, j) = pixels(width-1, j)*2+pixels_tmp((width-1)/2, j);
			}

			pixels_tmp = Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>(pixels);

			#pragma omp parallel for
			for(int i = 0; i < width; ++i)
			{
				for(int j = 0; j < height-1; j += 2)
				{
					pixels(i, j) = pixels_tmp(i, j/2);
				}
				for(int j = 1, j_2 = 0; j < height-1; j += 2, ++j_2)
				{
					float impair = pixels_tmp(i, width/2+j_2)*2;
					float pair_1 = pixels_tmp(i, j_2);
					float pair_2 = pixels_tmp(i, j_2+1);

					impair += (pair_1+pair_2)/2.f;
					pixels(i, j) = impair;
				}

				//Traitement spécifique pour la dernière ligne (différence avec le pair situé au dessus)
				pixels(i, height-1) = pixels(i, height-1)*2+pixels_tmp(i, (height-1)/2);
			}

			regenerateMap(mapOrigin, mapGenerated);
			project2DImageTo3DSpace(mapOrigin, mapGenerated);

			return true;
		}
	}
	return false;
}

bool Surface_DepthMapRendering_Plugin::savePointCloud(const QString& mapOrigin,
													  const QString& mapGenerated,
													  const QString& directory,
													  const int criteria,
													  const float radius)
{
	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(!directory.isEmpty() && mh_generated)
	{
		PFP2::MAP* generated_map = mh_generated->getMap();

		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("position");
		if(!position.isValid())
		{
			CGoGNerr << "position attribute is not valid" << CGoGNendl;
			return false;
		}
		VertexAttribute<PFP2::VEC3, PFP2::MAP> normal = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("normal");
		if(!position.isValid())
		{
			CGoGNerr << "normal attribute is not valid" << CGoGNendl;
			return false;
		}

		QString filename(directory);
		filename += "/" + mapOrigin + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += "PointClouds/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += QString::number(m_fbo->getWidth()) + "x" + QString::number(m_fbo->getHeight()) + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += mapGenerated;

		if(m_correspondance_done)
		{
			filename += "-Without-radius_"+QString::number(radius);

			switch(criteria)
			{
				case DENSITY:
					filename += "-Density.ply";
					break;
				case VISIBILITY:
					filename += "-Visibility.ply";
					break;
				default:
					break;
			}
		}
		else
		{
			filename += "-With.ply";
		}

		std::vector<VertexAttribute<typename PFP2::VEC3, typename PFP2::MAP>*> attribs;
		attribs.push_back(&position);
		attribs.push_back(&normal);

		return Algo::Surface::Export::exportPLYnew<PFP2>(*generated_map, attribs, filename.toStdString().c_str(), false);
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

//		pixels = (pixels.array()+1)/2.;

		int width = m_fbo->getWidth(), height = m_fbo->getHeight();

		QString filename(directory);
		filename += "/" + mapOrigin + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += "DepthMaps/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += QString::number(width) + "x" + QString::number(height) + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += mapGenerated;

		std::ofstream out;
		out.open(filename.toStdString() + "-originalDepthMap.dat", std::ios::out);
		if(!out.good())
		{
			CGoGNerr << "Unable to open file" << CGoGNendl;
			return false;
		}

		out << pixels;

//		for(int j = height-1; j >= 0; --j)
//		{
//			for(int i = 0; i < width; ++i)
//			{
//				out << pixels(i, j) << " " << std::flush;
//			}
//			out << std::endl;
//		}

		out.close();

		out.open(filename.toStdString() + "-MVPMatrix.dat", std::ios::out);
		if(!out.good())
		{
			CGoGNerr << "Unable to open file" << CGoGNendl;
			return false;
		}

		GLdouble mvp_matrix[16];
		camera->getModelViewProjectionMatrix(mvp_matrix);

		PFP2::MATRIX44 model_view_projection_matrix;

		for(int i = 0; i < 4; ++i)
		{
			for(int j = 0; j < 4; ++j)
			{
				model_view_projection_matrix(i, j) = mvp_matrix[i+j*4];
			}
		}

		out << model_view_projection_matrix;

		out.close();

		return true;
	}

	return false;
}

bool Surface_DepthMapRendering_Plugin::saveModifiedDepthMap(const QString& mapOrigin,
															const QString& mapGenerated,
															const QString& directory,
															const int criteria,
															const float radius)
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

		int width = m_fbo->getWidth(), height = m_fbo->getHeight();

		MapParameters& mapParams = m_mapParameterSet[mh_origin];
		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> pixels = mapParams.depthImageSet[mapGenerated];

		pixels = (pixels.array()+1)/2.;

		QString filename(directory);
		filename += "/" + mapOrigin + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += "DepthMaps/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += QString::number(width) + "x" + QString::number(height) + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += mapGenerated + "-modifiedDepthMap-radius_"+QString::number(radius);

		switch(criteria)
		{
			case DENSITY:
				filename += "-Density";
				break;
			case VISIBILITY:
				filename += "-Visibility";
				break;
			default:
				break;
		}

		QString filename2(filename);

		filename += ".dat";
		filename2 += "-Mask.dat";

		std::ofstream out;
		out.open(filename.toStdString(), std::ios::out);
		if(!out.good())
		{
			CGoGNerr << "Unable to open file" << CGoGNendl;
			return false;
		}

		std::ofstream out2;
		out2.open(filename2.toStdString(), std::ios::out);
		if(!out.good())
		{
			CGoGNerr << "Unable to open file" << CGoGNendl;
			return false;
		}

		Eigen::Matrix<GLint, Eigen::Dynamic, Eigen::Dynamic> mask_pixels;
		mask_pixels.setZero(pixels.rows(), pixels.cols());

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			int x = imageCoordinates[d].getXCoordinate(), y = imageCoordinates[d].getYCoordinate();

			if(fabs(1.f-pixels(x, y)) > FLT_EPSILON)
			{
				mask_pixels(x, y) = 1;
			}
		}

		for(int j = height-1; j >= 0; --j)
		{
			for(int i = 0; i < width; ++i)
			{
				if(mask_pixels(i, j)==1)
				{
					out << pixels(i, j) << " " << std::flush;
				}
				else
				{
					out << 1.f << " " << std::flush;
				}
				out2 << mask_pixels(i,j) << " " << std::flush;
			}
			out << std::endl;
			out2 << std::endl;
		}

		out.close();
		out2.close();

		return true;
	}

	return false;
}

bool Surface_DepthMapRendering_Plugin::saveMergedPointCloud(const QString& mapOrigin,
															const QStringList& mapNames,
															const QString& directory,
															const int criteria,
															const float radius)
{
	if(!mapOrigin.isEmpty() && !mapNames.empty() && !directory.isEmpty())
	{
		MapHandlerGen* mh_origin = m_schnapps->getMap(mapOrigin);
		if(m_mapParameterSet.contains(mh_origin))
		{
			MapParameters& mapParams = m_mapParameterSet[mh_origin];

			int level = *mapParams.decompositionLevelSet.begin();

			std::vector<PFP2::MAP*> maps;
			maps.reserve(mapNames.size());
			std::vector<std::vector<VertexAttribute<PFP2::VEC3, PFP2::MAP>>> attributes;
			for(int i = 0; i < mapNames.size(); ++i)
			{
				MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapNames[i]));
				if(mh_map)
				{
					std::vector<VertexAttribute<PFP2::VEC3, PFP2::MAP>> attribs;
					PFP2::MAP* map = mh_map->getMap();
					maps.push_back(map);

					VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
					VertexAttribute<PFP2::VEC3, PFP2::MAP> normal = mh_map->getAttribute<PFP2::VEC3, VERTEX>("normal");
					attribs.push_back(position);
					attribs.push_back(normal);
					attributes.push_back(attribs);
				}
				else
				{
					return false;
				}
			}

			int width = m_fbo->getWidth(), height = m_fbo->getHeight();

			QString filename(directory);
			filename += "/" + mapOrigin + "/";
			mkdir(filename.toStdString().c_str(), 0777);

			filename += "PointClouds/";
			mkdir(filename.toStdString().c_str(), 0777);

			filename += QString::number(width) + "x" + QString::number(height) + "/";
			mkdir(filename.toStdString().c_str(), 0777);

			if(m_correspondance_done)
			{
				filename += mapOrigin + "-" + QString::number(width) + "x" + QString::number(height) + "-Merged-Without-Level-" + QString::number(level)
						+"-radius_"+QString::number(radius);

				switch(criteria)
				{
					case DENSITY:
						filename += "-Density.ply";
						break;
					case VISIBILITY:
						filename += "-Visibility.ply";
						break;
					default:
						break;
				}
			}
			else
			{
				filename += mapOrigin + "-" + QString::number(width) + "x" + QString::number(height) + "-Merged-With-Level-" + QString::number(level) + ".ply";
			}

			return Algo::Surface::Export::exportPLYVertMaps<PFP2>(maps, attributes, filename.toStdString().c_str(), false);
		}
	}

	return false;
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

		const int width = m_fbo->getWidth(), height = m_fbo->getHeight();

		Eigen::Matrix<GLfloat, Eigen::Dynamic, 3> matrix(width*height, 3);

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			matrix.row(imageCoordinates[d].getXCoordinate()*height+imageCoordinates[d].getYCoordinate())
					= Eigen::Vector3f(position[d][0], position[d][1], position[d][2]);
		}

		for(int i = 0; i < matrix.rows()-height; ++i)
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

		mh_generated->notifyAttributeModification(normal, false);

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

		qglviewer::Vec vd = -(camera->viewDirection());
		PFP2::VEC3 view_direction(vd.x, vd.y, vd.z);

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>& depthImage = mapParams.depthImageSet[mapGenerated];

		VertexAttribute<PFP2::VEC3, PFP2::MAP> normal = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("normal");
		VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_generated->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
		VertexAttribute<float, PFP2::MAP> visibilityConfidence = mh_generated->addAttribute<float, VERTEX>("VisibilityConfidence");
		VertexAttribute<float, PFP2::MAP> label = mh_generated->addAttribute<float, VERTEX>("Label");

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> gradientMagnitude;
		gradientMagnitude.setZero(depthImage.rows(), depthImage.cols());

		//Calcul de l'intensité du gradient (sqrt(grad_x^2+grad_y^2))
//		#pragma omp parallel for
//		for(int i = 1; i < depthImage.rows()-1; ++i)
//		{
//			for(int j = 1; j < depthImage.cols()-1; ++j)
//			{
//				float grad_x = (depthImage(i+1,j)-depthImage(i-1,j))/2.f;
//				float grad_y = (depthImage(i,j+1)-depthImage(i,j-1))/2.f;
//				gradientMagnitude(i, j) = sqrt(grad_x*grad_x+grad_y*grad_y);
//			}
//		}

//		const float mean = gradientMagnitude.mean();

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			int x = imageCoordinates[d].getXCoordinate(), y = imageCoordinates[d].getYCoordinate();
			float color = depthImage(x, y);
			if(fabs(1-color) > FLT_EPSILON/* && gradientMagnitude(x, y) < mean+mean/2.f*/)
			{
//				PFP2::VEC3 u = (position_camera-position[d]).normalized();
				visibilityConfidence[d] = view_direction*normal[d];
				if(visibilityConfidence[d] != visibilityConfidence[d])
				{	//visibilityConfidence[d] is NaN
					visibilityConfidence[d] = 0.f;
				}
				label[d] = d.label();
			}
			else
			{
				//Suppression du point de la carte de profondeur
				depthImage(x, y) = 1.f;
				label[d] = Dart::nil().label();
			}
		}

		CGoGNout << ".. fait en " << chrono.elapsed() << " ms" << CGoGNendl;

		deleteBackground(mapOrigin, mapGenerated);

		mh_generated->notifyAttributeModification(visibilityConfidence, false);
		mh_generated->notifyAttributeModification(label, false);
	}
}

void Surface_DepthMapRendering_Plugin::densityEstimation(const QString& mapOrigin, const QString& mapGenerated, const float radius)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin))
	{
		CGoGNout << "Estimation de la densité des points de la carte " << mapGenerated.toStdString() << " .." << CGoGNflush;
		Utils::Chrono chrono;
		chrono.start();

		MapParameters& mapParams = m_mapParameterSet[mh_origin];
		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>& depthImage = mapParams.depthImageSet[mapGenerated];

		PFP2::MAP* generated_map = mh_generated->getMap();
		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("position");
		VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_generated->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
		VertexAttribute<float, PFP2::MAP> samplingDensity = mh_generated->getAttribute<float, VERTEX>("SamplingDensity");
		if(!samplingDensity.isValid())
		{
			samplingDensity = mh_generated->addAttribute<float, VERTEX>("SamplingDensity");
		}

		Eigen::Matrix<Dart, Eigen::Dynamic, Eigen::Dynamic> imageDarts;
		imageDarts.setZero(depthImage.rows(), depthImage.cols());

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			int x = imageCoordinates[d].getXCoordinate(), y = imageCoordinates[d].getYCoordinate();
			imageDarts(x, y) = d;
		}

		//Rayon normalisé par rapport à la diagonale de la boîte englobante du nuage de points original
		const float NORM_RADIUS_SPHERE = radius*mh_origin->getBBdiagSize();
		//Au carré pour n'avoir qu'à calculer la norme au carré
		const float RADIUS2_SPHERE = NORM_RADIUS_SPHERE*NORM_RADIUS_SPHERE;
		const float VOLUME_SPHERE = (4./3.)*M_PI*(NORM_RADIUS_SPHERE*NORM_RADIUS_SPHERE*NORM_RADIUS_SPHERE);

//		std::cout << std::endl << NORM_RADIUS_SPHERE << std::endl;

		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			int x = imageCoordinates[d].getXCoordinate(), y = imageCoordinates[d].getYCoordinate();
			int count = 1;	//Toujours 1 au début (parce qu'il y a au moins un point appartenant à la sphère : lui-même

			bool stop_search = false;
			int neighborhood = 1;
			while(!stop_search)
			{
				stop_search = true;

				bool do_min_i = true, do_max_i = true;

				int min_i = x-neighborhood;
				if(min_i<0)
				{
					do_min_i = false;
					min_i = 0;
				}

				int max_i = x+neighborhood;
				if(max_i>=depthImage.cols())
				{
					do_max_i = false;
					max_i = depthImage.cols()-1;
				}

				bool do_min_j = true, do_max_j = true;

				int min_j = y-neighborhood;
				if(min_j<0)
				{
					do_min_j = false;
					min_j = 0;
				}

				int max_j = y-neighborhood;
				if(max_j>=depthImage.rows())
				{
					do_max_j = false;
					max_j = depthImage.rows()-1;
				}

//				#pragma omp parallel for shared(count)
				for(int i = min_i; i <= max_i; ++i)
				{
					if(do_min_j && fabs(1.f - depthImage(i, min_j)) > FLT_EPSILON)
					{
						if((position[d]-position[imageDarts(i, min_j)]).norm2() < RADIUS2_SPHERE)	//Comparaison avec la norme du vecteur au carré (plus rapide)
						{
							//Si le point est contenu dans la sphere de rayon NORM_RADIUS_SPHERE
							++count;
							stop_search = false;
						}
					}

					if(do_max_j && fabs(1.f - depthImage(i, max_j)) > FLT_EPSILON)
					{
						if((position[d]-position[imageDarts(i, max_j)]).norm2() < RADIUS2_SPHERE)	//Comparaison avec la norme du vecteur au carré (plus rapide)
						{
							//Si le point est contenu dans la sphere de rayon NORM_RADIUS_SPHERE
							++count;
							stop_search = false;
						}
					}
				}

				++min_j;
				--max_j;

//				#pragma omp parallel for shared(count)
				for(int j = min_j; j <= max_j; ++j)
				{
					if(do_min_i && fabs(1.f - depthImage(min_i, j)) > FLT_EPSILON)
					{
						if((position[d]-position[imageDarts(min_i, j)]).norm2() < RADIUS2_SPHERE)	//Comparaison avec la norme du vecteur au carré (plus rapide)
						{
							//Si le point est contenu dans la sphere de rayon NORM_RADIUS_SPHERE
							++count;
							stop_search = false;
						}
					}

					if(do_max_i && fabs(1.f - depthImage(max_i, j)) > FLT_EPSILON)
					{
						if((position[d]-position[imageDarts(max_i, j)]).norm2() < RADIUS2_SPHERE)	//Comparaison avec la norme du vecteur au carré (plus rapide)
						{
							//Si le point est contenu dans la sphere de rayon NORM_RADIUS_SPHERE
							++count;
							stop_search = false;
						}
					}
				}
				++neighborhood;
			}
			samplingDensity[d] = count/VOLUME_SPHERE;
		}

		CGoGNout << ".. fait en " << chrono.elapsed() << " ms" << CGoGNendl;
	}
}

void Surface_DepthMapRendering_Plugin::findCorrespondingPoints(const QString& mapOrigin, const QString& mapGenerated, int criteria)
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
		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>& depthImage = mapParams.depthImageSet[mapGenerated];
		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>& newDepthImage = mapParams.newDepthImageSet[mapGenerated];
		newDepthImage = Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>(depthImage);

		PFP2::MAP* generated_map = mh_generated->getMap();
		VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("position");
		VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_generated->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
//		VertexAttribute<NoTypeNameAttribute<std::vector<PointCorrespondance>>, PFP2::MAP> correspondingPointsAttribute
//				= mh_generated->addAttribute<NoTypeNameAttribute<std::vector<PointCorrespondance>>, VERTEX>("CorrespondingPoints");
		const int width = m_fbo->getWidth(), height = m_fbo->getHeight();


		VertexAttribute<float, PFP2::MAP> criteriaAttribute;

		switch(criteria)
		{
			case DENSITY:
				criteriaAttribute = mh_generated->getAttribute<float, VERTEX>("SamplingDensity");
				break;
			case VISIBILITY:
				criteriaAttribute = mh_generated->getAttribute<float, VERTEX>("VisibilityConfidence");
				break;
			default:
				break;
		}

		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> currentDepthImage;
		currentDepthImage.setOnes(width, height);
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> labelValues;
		labelValues.setZero(width, height);
		Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> maxConfidenceValues;
		maxConfidenceValues.setZero(width, height);
		Eigen::Matrix<Dart, Eigen::Dynamic, Eigen::Dynamic> imageDarts;
		imageDarts.setZero(width, height);

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			int x = imageCoordinates[d].getXCoordinate(), y = imageCoordinates[d].getYCoordinate();
			imageDarts(x, y) = d;
		}

		//Change current rendering camera, and reposition it correctly
		qglviewer::Vec camera_position = camera->position();
		m_schnapps->getSelectedView()->setCurrentCamera(camera);
		camera->setPosition(camera_position);

		float z_near = camera->zNear(), z_far = camera->zFar();

		/*
		 * z_w = valeur de profondeur de l'image ; dans l'intervalle [0;1]
		 * n = distance non signée du zNear à la caméra
		 *
		 * f = distance non signée du zFar à la caméra
		 * p = (f*n)/(z_w*(f-n)-f) => Fonction qui calcule la distance d'un point à la caméra, selon la valeur de profondeur de la carte de profondeur
		 * - ((f*n) * (f-n)) / ((z_w*(f-n)-f)*(z_w*(f-n)-f)) => Dérivée 1ère de la fonction précédente
		 * z_w = (f*(p+n))/(p*(f-n)) => Fonction inverse qui pour une distance d'un point à la caméra, donne sa valeur de profondeur
		 */

		const float threshold = 1/64.f;

//		std::vector<Dart> tmp_vector;
//		tmp_vector.resize(11);

//		std::vector<std::vector<Dart>> correspondingPoints;
//		correspondingPoints.resize(width*height, tmp_vector);

//		std::vector<MapHandlerGen*> vec_maps;
//		vec_maps.reserve(11);

		for(QHash<QString, MapHandlerGen*>::iterator it = mapParams.projectedMapSet.begin(); it != mapParams.projectedMapSet.end(); ++it)
		{
			if(it.value() != mhg_generated)
			{
				MapHandler<PFP2>* mh_current = static_cast<MapHandler<PFP2>*>(it.value());

				m_shaderScalarFieldReal->setAttributePosition(mh_current->getVBO("position"));
				m_shaderScalarFieldReal->setAttributeScalar(mh_current->getVBO("Label"));

				VertexAttribute<float, PFP2::MAP> criteriaAttributeCurrent;

				switch(criteria)
				{
					case DENSITY:
						criteriaAttributeCurrent = mh_current->getAttribute<float, VERTEX>("SamplingDensity");
						break;
					case VISIBILITY:
						criteriaAttributeCurrent = mh_current->getAttribute<float, VERTEX>("VisibilityConfidence");
						break;
					default:
						break;
				}

//				VertexAttribute<PFP2::VEC3, PFP2::MAP> positionCurrent = mh_current->getAttribute<PFP2::VEC3, VERTEX>("position");

				m_fbo->bind();
				glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );	//To clean the color and depth textures

				mh_current->draw(m_shaderScalarFieldReal, CGoGN::Algo::Render::GL2::POINTS);

				glEnable(GL_TEXTURE_2D);

				glBindTexture(GL_TEXTURE_2D, *m_fbo->getDepthTexId());
				glGetTexImage(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, GL_FLOAT, currentDepthImage.data());

				glBindTexture(GL_TEXTURE_2D, *m_fbo->getColorTexId(0));
				glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, labelValues.data());

				glDisable(GL_TEXTURE_2D);
				m_fbo->unbind();

//				Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic> copyImage(currentDepthImage);

//				currentDepthImage = currentDepthImage.array()*2-1;    //Set range to [-1;1]

				currentDepthImage = (depthImage.array()-currentDepthImage.array()).abs();

//				currentDepthImage = (z_far*z_near)/(((currentDepthImage.array()+1)/2)*(z_far-z_near)+z_far);

//				int count = 0, count_0 = 0;

//				for(int i = 0; i < currentDepthImage.rows(); ++i)
//				{
//					for(int j = 0; j < currentDepthImage.cols(); ++j)
//					{

//						if(fabs(1-depthImage(i, j)) < FLT_EPSILON)
//						{
//							++count_0;
//						}
//						else if(currentDepthImage(i, j) < threshold && fabs(1-depthImage(i, j)) > FLT_EPSILON)
//						{
//							++count;
//						}
//					}
//				}

//				CGoGNout << count << " point(s) similaire(s) sur " << depthImage.rows()*depthImage.cols()-count_0 << " point(s) au total." << CGoGNendl;

//				QImage image(copyImage.rows(), copyImage.cols(), QImage::Format_RGB32);
//				for(int i = 0; i < currentDepthImage.rows(); ++i)
//				{
//					for(int j = 0; j < currentDepthImage.cols(); ++j)
//					{
//						if(currentDepthImage(i, j) < threshold)
//						{
//							image.setPixel(i, copyImage.rows()-j-1, qRgb((1-copyImage(i,j))*255, (1-copyImage(i,j))*255, (1-copyImage(i,j))*255));
//						}
//						else
//						{
//							image.setPixel(i, copyImage.rows()-j-1, qRgb(0, 0, 0));
//						}
//					}
//				}

//				image.save("/home/blettere/Projets/Test/Images/"+mapGenerated+"-"+mh_current->getName()+".png");

//				unsigned int size = vec_maps.size();

				#pragma omp parallel for
				for(int i = 0; i < currentDepthImage.rows(); ++i)
				{
					for(int j = 0; j < currentDepthImage.cols(); ++j)
					{
						//Si un point s'est projeté sur ce pixel
						Dart d = Dart::create(labelValues(i, j));
						if(d != Dart::nil() && d.label() < 4000000000)
						{
							//Dummy test to avoid problems

//							correspondingPoints[i+j*width][size] = Dart::nil();

//							float d_o = (z_far*z_near)/(((depthImage(i, j)+1)/2)*(z_far-z_near)+z_far);

//							CGoGNout << fabs(d_o-currentDepthImage(i, j)) << CGoGNendl;

//							if(fabs(1-depthImage(i, j)) > FLT_EPSILON && fabs(d_o-currentDepthImage(i, j)) < threshold)
							if(fabs(1-depthImage(i, j)) > FLT_EPSILON && currentDepthImage(i, j) < threshold)
							{
//								correspondingPoints[i+j*width][size] = d;
								maxConfidenceValues(i, j) = std::max(maxConfidenceValues(i, j), criteriaAttributeCurrent[d]);
							}
						}
					}
				}
//				vec_maps.push_back(it.value());
			}
		}

		m_schnapps->getSelectedView()->setCurrentCamera("camera_0");

		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			int x = imageCoordinates[d].getXCoordinate(), y = imageCoordinates[d].getYCoordinate();
//			correspondingPointsAttribute[d].reserve(11);
//			for(unsigned int i = 0; i < vec_maps.size(); ++i)
//			{
//				if(vec_maps[i] != mhg_generated && correspondingPoints[x+y*width][i] != Dart::nil())
//				{
//					PointCorrespondance tmp;
//					tmp.map = vec_maps[i];
//					tmp.vertex = correspondingPoints[x+y*width][i];
//					correspondingPointsAttribute[d].push_back(tmp);
//				}
//			}
			if(criteriaAttribute[d] < maxConfidenceValues(x, y))
			{
				//Suppression du point dans la nouvelle carte de profondeur (pour ne pas influer les comparaisons futures)
				newDepthImage(x, y) = 1.f;
//				depthImage(x, y) = 1.f;
			}
		}
//		deleteBackground(mapOrigin, mapGenerated);
		mh_generated->notifyAttributeModification(criteriaAttribute, false);
	}
}

void Surface_DepthMapRendering_Plugin::updateDepthImages(const QString& mapOrigin)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	if(mh_origin && m_mapParameterSet.contains(mh_origin))
	{
		MapParameters& mapParams = m_mapParameterSet[mh_origin];
		for(QHash<QString, Eigen::Matrix<GLfloat, Eigen::Dynamic, Eigen::Dynamic>>::iterator it = mapParams.newDepthImageSet.begin(); it != mapParams.newDepthImageSet.end(); ++it)
		{
			const QString mapGenerated = it.key();
			mapParams.depthImageSet[mapGenerated] = it.value();
			it.value().resize(0, 0);

			deleteBackground(mapOrigin, mapGenerated);
		}
	}

	m_correspondance_done = true;
}

void Surface_DepthMapRendering_Plugin::regenerateMap(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin))
	{
		PFP2::MAP* generated_map = mh_generated->getMap();

		generated_map->clear(false);

		VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinates = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
		if(!planeCoordinates.isValid())
		{
			planeCoordinates = mh_generated->addAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
		}
		VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinates = mh_generated->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");

		MapParameters& mapParams = m_mapParameterSet[mh_origin];

		const int level = mapParams.decompositionLevelSet[mapGenerated];
		const int img_width = m_fbo->getWidth(), img_height = m_fbo->getHeight();

		const int l_p = pow(2, level);
		const int width = img_width/l_p;
		const int height = img_height/l_p;

		Algo::Surface::Tilings::Square::Grid<PFP2> grid(*generated_map, width-1, height-1);	//Create a grid of width*height vertices
		grid.embedIntoGrid(planeCoordinates, 2, 2);

		mh_generated->updateBB(planeCoordinates);

		qglviewer::Vec bb_min = mh_generated->getBBmin(), bb_max = mh_generated->getBBmax();

		const float step_width = (bb_max.x-bb_min.x)/img_width, step_height = (bb_max.y-bb_min.y)/img_height;

		const float shift_width = (l_p-1)*step_width, shift_height = (l_p-1)*step_height;

		grid.embedIntoGrid(planeCoordinates, 2-shift_width, 2-shift_height);

		const float shift_width_2 = shift_width/2.f, shift_height_2 = shift_height/2.f;

		std::vector<Dart>& vDarts = grid.getVertexDarts();

		for(int j = 0; j < height; ++j)
		{
			for(int i = 0; i < width; ++i)
			{
				Dart d = vDarts[j*width+i];
				planeCoordinates[d][0] -= shift_width_2;
				planeCoordinates[d][1] -= shift_height_2;

				imageCoordinates[d].setCoordinates(i, j);
			}
		}

		mh_generated->notifyConnectivityModification(false);
		mh_generated->notifyAttributeModification(planeCoordinates, false);
		mh_generated->notifyAttributeModification(imageCoordinates, false);
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
					generated_map->deleteFace(d);
					stop = true;
				}
			}
		}

		mh_generated->notifyConnectivityModification(false);
		if(position.isValid())
		{
			mh_generated->updateVBO(position);
			mh_generated->updateBB(position);
		}

		mh_generated->notifyConnectivityModification(false);
		mh_generated->notifyAttributeModification(position, false);
		mh_generated->notifyAttributeModification(imageCoordinates, false);
	}
}

void Surface_DepthMapRendering_Plugin::removeUselessAttributes(const QString &mapGenerated)
{
	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_generated)
	{
		PFP2::MAP* generated_map = mh_generated->getMap();

		VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinates = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
		VertexAttribute<PFP2::VEC3, PFP2::MAP> normal = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("normal");
		VertexAttribute<float, PFP2::MAP> visibilityConfidence = mh_generated->getAttribute<float, VERTEX>("VisibilityConfidence");
		VertexAttribute<float, PFP2::MAP> label = mh_generated->getAttribute<float, VERTEX>("Label");

		generated_map->removeAttribute(planeCoordinates);
		generated_map->removeAttribute(normal);
		generated_map->removeAttribute(visibilityConfidence);
		generated_map->removeAttribute(label);
	}
}

void Surface_DepthMapRendering_Plugin::verifyDepthMaps(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);

	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin))
	{
		//Iterate through camera set, verifying the alignment of the camera plane (viewing direction) with the current camera plane (viewing direction)
	}
}

#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_Plugin, Surface_DepthMapRendering_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_PluginD, Surface_DepthMapRendering_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
