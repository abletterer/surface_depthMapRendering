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

	m_shader = new CGoGN::Utils::ShaderSimpleColor();

	registerShader(m_shader);

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
	if(m_shader)
	{
		delete m_shader;
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
		MapParameters& mapParam = i.value();
		if(mapParam.positionVBO == vbo)
		{
			mapParam.positionVBO = NULL;
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

void Surface_DepthMapRendering_Plugin::createCameras(const QString& mapName)
{
	MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
	MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);

	if(mh_map && m_mapParameterSet.contains(mhg_map))
	{
		QString baseName("DepthCamera-");

		MapParameters& mapParam = m_mapParameterSet[mhg_map];

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

		for(int i = 0; i < 12; ++i)
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

			mapParam.depthCameraSet[generatedName] = camera;
		}
	}
}

void Surface_DepthMapRendering_Plugin::render(const QString& mapName, const QString& directory)
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

		MapParameters& mapParam = m_mapParameterSet[mhg_map];

		int width = m_depthFBO->getWidth(), height = m_depthFBO->getHeight();
		m_shader->setAttributePosition(mapParam.positionVBO);

		std::vector<GLfloat> pixels;
		pixels.resize(width*height);

		Utils::Chrono chrono;
		chrono.start();

		std::vector<QString> mapNames;
		mapNames.reserve(mapParam.depthCameraSet.size());

		for(QHash<QString, Camera*>::iterator it = mapParam.depthCameraSet.begin(); it != mapParam.depthCameraSet.end(); ++it)
		{
			Camera* camera = it.value();
			QString cameraName(camera->getName());

			QString generatedName(mapName);
			generatedName += "-" + cameraName;

			mapNames.push_back(generatedName);

			m_schnapps->getSelectedView()->setCurrentCamera(camera);

			m_depthFBO->bind();
			glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );	//To clean the color and depth textures
			mh_map->draw(m_shader, CGoGN::Algo::Render::GL2::TRIANGLES);	//Render the map into the FrameBufferObject

			glBindTexture(GL_TEXTURE_2D, *m_depthFBO->getDepthTexId());

			//Read pixels of the generated texture and store them in an array
			glGetTexImage(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, GL_FLOAT, pixels.data());
			m_depthFBO->unbind();

			m_schnapps->getSelectedView()->setCurrentCamera("camera_0");

			MapHandlerGen* mhg_generated = m_schnapps->addMap(generatedName, 2);
			mapParam.projectedMapSet[generatedName] = mhg_generated;

			mapParam.decompositionLevel[generatedName] = 0;

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
					pixels[i+width*j] = 2*pixels[i+width*j]-1;
				}
			}

			mapParam.depthImageSet[generatedName] = pixels;

			mh_generated->notifyConnectivityModification(false);
			mh_generated->notifyAttributeModification(planeCoordinatesGenerated, false);
			mh_generated->notifyAttributeModification(imageCoordinatesGenerated, false);
			project2DImageTo3DSpace(mapName, generatedName);
		}

		CGoGNout << "Temps d'échantillonnage : " << chrono.elapsed() << " ms " << CGoGNflush;
		CGoGNout << "pour " << mapParam.depthCameraSet.size() << " vue(s) différente(s) " << CGoGNflush;
		CGoGNout << "de taille " << width << "x" << height << CGoGNflush;
		CGoGNout << " sur un objet composé de " << mh_map->getMap()->getNbCells(VERTEX) << " point(s)" << CGoGNendl;

		saveMergedPointCloud(mapName, mapNames, directory);

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
		MapParameters& mapParam = m_mapParameterSet[mhg_origin];
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

		Camera* camera = mapParam.depthCameraSet[mapGenerated];
		std::vector<GLfloat> pixels = mapParam.depthImageSet[mapGenerated];

		TraversorF<PFP2::MAP> trav_face_map(*generated_map);
		Dart next;
		bool stop = false;
		for(Dart d = trav_face_map.begin(); d != trav_face_map.end(); d = next)
		{

			next = trav_face_map.next();
			stop = false;
			Traversor2FV<PFP2::MAP> trav_vert_face_map(*generated_map, d);
			for(Dart dd = trav_vert_face_map.begin(); !stop && dd != trav_vert_face_map.end(); dd = trav_vert_face_map.next())
			{
				float color = pixels[imageCoordinatesGenerated[dd].getXCoordinate()+m_depthFBO->getWidth()*imageCoordinatesGenerated[dd].getYCoordinate()];
				if(fabs(1-color)<FLT_EPSILON)
				{
					generated_map->deleteFace(d);
					stop = true;
				}
			}
		}

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
			float color = pixels[imageCoordinatesGenerated[d].getXCoordinate()+m_depthFBO->getWidth()*imageCoordinatesGenerated[d].getYCoordinate()];

			PFP2::VEC4 pos = PFP2::VEC4(planeCoordinatesGenerated[d][0], planeCoordinatesGenerated[d][1], color, 1.f);

			pos = model_view_projection_matrix_inv*pos;

			positionGenerated[d] = PFP2::VEC3(pos[0]/pos[3], pos[1]/pos[3], pos[2]/pos[3]);
		}

		mh_generated->updateBB(positionGenerated);
		mh_generated->notifyAttributeModification(positionGenerated, false);
	}
}

bool Surface_DepthMapRendering_Plugin::moveDownDecomposition(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapOrigin));
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(m_schnapps->getMap(mapGenerated));

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mh_origin))
	{
		MapParameters& mapParam = m_mapParameterSet[mh_origin];
		PFP2::MAP* generated_map = mh_generated->getMap();

		VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinatesGenerated = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
		if(!planeCoordinatesGenerated.isValid())
		{
			CGoGNerr << "PlaneCoordinates attribute is not valid" << CGoGNendl;
			return false;
		}

		VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinatesGenerated = mh_generated->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
		if(!imageCoordinatesGenerated.isValid())
		{
			CGoGNerr << "ImageCoordinates attribute is not valid" << CGoGNendl;
			return false;
		}

		int img_width = m_depthFBO->getWidth();
		int img_height = m_depthFBO->getHeight();

		int& level = mapParam.decompositionLevel[mapGenerated];

		int l_p = pow(2, level);

		int current_width = img_width/l_p;
		int current_height = img_height/l_p;

		if(current_width < 8 || current_height < 8)
		{
			return false;
		}

		generated_map->clear(false);

		std::vector<GLfloat>& matrix = mapParam.depthImageSet[mapGenerated];

		std::vector<GLfloat> matrix2 = std::vector<GLfloat>(matrix);

		for(int i = 0; i < current_width; ++i)
		{
			for(int j = 0; j < current_height; ++j)
			{
				if(i%2==0)
				{
					matrix[i/2+img_width*j] = matrix2[i+img_width*j];
				}
				else
				{
					int left = matrix2[i-1+img_width*j];
					int result = matrix2[i+img_width*j];
					if(i != current_width-1)
					{
						int right = matrix2[i+1+img_width*j];
						result -= (left+right)/2.f;
						matrix[current_width/2+i/2+img_width*j] = result;
					}
					else
					{
						result -= left;
						matrix[current_width/2+i/2+img_width*j] = result;
					}
				}
			}
		}

		matrix2 = std::vector<GLfloat>(matrix);

		for(int i = 0; i < current_width; ++i)
		{
			for(int j = 0; j < current_height; ++j)
			{
				if(j%2==0)
				{
					matrix[i+img_width*(j/2)] = matrix2[i+img_width*j];
				}
				else
				{
					int top = matrix2[i+img_width*(j-1)];
					int result = matrix2[i+img_width*j];
					if(j != current_height-1)
					{
						int bottom = matrix2[i+img_width*(j+1)];
						result -= (top+bottom)/2.f;
						matrix[i+img_width*(current_height/2+j/2)] = result;
					}
					else
					{
						result -= top;
						matrix[i+img_width*(current_height/2+j/2)] = result;
					}
				}
			}
		}

		++level;

		current_width /= 2;
		current_height /= 2;

		Algo::Surface::Tilings::Square::Grid<PFP2> grid(*generated_map, current_width-1, current_height-1);
		grid.embedIntoGrid(planeCoordinatesGenerated, img_width-1, img_height-1);

		mh_generated->updateBB(planeCoordinatesGenerated);

		qglviewer::Vec bb_min = mh_generated->getBBmin();
		qglviewer::Vec bb_max = mh_generated->getBBmax();

		float width_step = (bb_max.x-bb_min.x)/img_width;
		float height_step = (bb_max.y-bb_min.y)/img_height;

		grid.embedIntoGrid(planeCoordinatesGenerated, img_width-l_p, img_height-l_p);

		PFP2::MATRIX44 transform_matrix;
		transform_matrix.identity();
		transform_matrix.setSubVectorV(0, 3, PFP2::VEC4(-(width_step*(l_p-1))/2.f, (height_step*(l_p-1))/2.f, 0., 1.));
		grid.transform(planeCoordinatesGenerated, transform_matrix);

		std::vector<Dart>& vDarts = grid.getVertexDarts();

		for(int i = 0; i < current_width; ++i)
		{
			for(int j = 0; j < current_height; ++j)
			{
				//Set plane coordinates in [-1;1]
				planeCoordinatesGenerated[vDarts[j*current_width+i]][0] /= (img_width-l_p)/2.f;
				planeCoordinatesGenerated[vDarts[j*current_width+i]][1] /= (img_height-l_p)/2.f;

				imageCoordinatesGenerated[vDarts[i+current_width*j]].setCoordinates(i, j);
			}
		}

		mh_generated->notifyAttributeModification(planeCoordinatesGenerated, false);
		mh_generated->notifyAttributeModification(imageCoordinatesGenerated, false);
		mh_generated->notifyConnectivityModification(false);

		project2DImageTo3DSpace(mapOrigin, mapGenerated);
		m_schnapps->getSelectedView()->updateGL();
	}
	else
	{
		return false;
	}
	return true;
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

		QString filename(directory);
		filename += "/" + mapOrigin + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += "DepthMaps/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += QString::number(m_depthFBO->getWidth()) + "x" + QString::number(m_depthFBO->getHeight()) + "/";
		mkdir(filename.toStdString().c_str(), 0777);

		filename += mapGenerated + ".dat";

		std::ofstream out;
		out.open(filename.toStdString(), std::ios::out);
		if(!out.good())
		{
			CGoGNerr << "Unable to open file" << CGoGNendl;
			return false;
		}

		std::vector<GLfloat> depthMap = mapParams.depthImageSet[mapGenerated];

		for(int j = m_depthFBO->getHeight()-1; j >= 0; --j)
		{
			for(int i = 0; i < m_depthFBO->getWidth(); ++i)
			{
				out << depthMap[i+m_depthFBO->getWidth()*j] << " " << std::flush;
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

#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_Plugin, Surface_DepthMapRendering_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_PluginD, Surface_DepthMapRendering_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
