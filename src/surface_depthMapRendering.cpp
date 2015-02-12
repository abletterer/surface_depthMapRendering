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

	m_depthShader = new CGoGN::Utils::ShaderDepth();

	m_depthFBO = new CGoGN::Utils::FBO(256, 256);
	m_depthFBO->createAttachDepthTexture();
	m_depthFBO->createAttachColorTexture(GL_RGBA);

	registerShader(m_depthShader);

	connect(m_depthMapRenderingAction, SIGNAL(triggered()), this, SLOT(openDepthMapRenderingDialog()));

	connect(m_depthMapRenderingDialog, SIGNAL(accepted()), this, SLOT(renderFromDialog()));
	connect(m_depthMapRenderingDialog->button_cancel, SIGNAL(clicked()), this, SLOT(closeDepthMapRenderingDialog()));
	connect(m_depthMapRenderingDialog->button_ok, SIGNAL(clicked()), this, SLOT(renderFromDialog()));

	connect(m_schnapps, SIGNAL(mapAdded(MapHandlerGen*)), this, SLOT(mapAdded(MapHandlerGen*)));
	connect(m_schnapps, SIGNAL(mapRemoved(MapHandlerGen*)), this, SLOT(mapRemoved(MapHandlerGen*)));

	foreach(MapHandlerGen* map, m_schnapps->getMapSet().values())
		mapAdded(map);

	return true;
}

void Surface_DepthMapRendering_Plugin::disable()
{
	if(m_depthShader)
	{
		delete m_depthShader;
	}
	if(m_depthFBO)
	{
		delete m_depthFBO;
	}

	disconnect(m_depthMapRenderingAction, SIGNAL(triggered()), this, SLOT(openDepthMapRenderingDialog()));

	disconnect(m_depthMapRenderingDialog, SIGNAL(accepted()), this, SLOT(renderFromDialog()));
	disconnect(m_depthMapRenderingDialog->button_cancel, SIGNAL(clicked()), this, SLOT(closeDepthMapRenderingDialog()));
	disconnect(m_depthMapRenderingDialog->button_ok, SIGNAL(clicked()), this, SLOT(renderFromDialog()));

	disconnect(m_schnapps, SIGNAL(mapAdded(MapHandlerGen*)), this, SLOT(mapAdded(MapHandlerGen*)));
	disconnect(m_schnapps, SIGNAL(mapRemoved(MapHandlerGen*)), this, SLOT(mapRemoved(MapHandlerGen*)));
}

void Surface_DepthMapRendering_Plugin::drawMap(View* view, MapHandlerGen* map)
{
	if(m_depthShader && m_mapParameterSet.contains(map))
	{
		float near = view->getCurrentCamera()->zNear();
		float far = view->getCurrentCamera()->zFar();

		m_depthShader->setZMin(near/far);
		m_depthShader->setZMax(1.);
		m_depthShader->setAttributePosition(m_mapParameterSet[map].positionVBO);

		map->draw(m_depthShader, CGoGN::Algo::Render::GL2::TRIANGLES);
	}
}

void Surface_DepthMapRendering_Plugin::openDepthMapRenderingDialog()
{
	m_depthMapRenderingDialog->show();
}

void Surface_DepthMapRendering_Plugin::closeDepthMapRenderingDialog()
{
	m_depthMapRenderingDialog->close();
}

void Surface_DepthMapRendering_Plugin::renderFromDialog()
{
	QList<QListWidgetItem*> currentItems = m_depthMapRenderingDialog->list_maps->selectedItems();
	if(!currentItems.empty())
	{
		render(currentItems[0]->text());
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

void Surface_DepthMapRendering_Plugin::createCameras(const QString& mapName, int number)
{
	MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
	MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);

	if(mh_map && m_mapParameterSet.contains(mhg_map))
	{
		QString baseName("DepthCamera-");

		MapParameters& mapParam = m_mapParameterSet[mhg_map];

		for(int i = 0; i < number; ++i)
		{
			QString cameraName(baseName);
			cameraName.append(QString::number(i));
			Camera* camera = m_schnapps->addCamera(cameraName);

			qglviewer::Vec bb_min = mh_map->getBBmin();
			qglviewer::Vec bb_max = mh_map->getBBmax();

			qglviewer::Vec center = (bb_min+bb_max)/2.f;

			camera->setSceneBoundingBox(bb_min,bb_max);
			camera->showEntireScene();

			qglviewer::Vec camera_position(camera->position());

			float radius = qAbs(camera_position.z - center.z);

			camera_position.x = center.x + radius*std::cos(M_PI/(number/2)*i);
			camera_position.y = center.y;
			camera_position.z = center.z + radius*std::sin(M_PI/(number/2)*i);

			camera->setPosition(camera_position);

			camera->lookAt(center);

			QString generatedName(mapName);
			generatedName.append("-");
			generatedName.append(cameraName);

			mapParam.depthCameraSet[generatedName] = camera;
		}
	}
}

void Surface_DepthMapRendering_Plugin::render(const QString& mapName, const QString& directory)
{
	if(m_depthShader)
	{
		MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
		MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);

		if(mh_map && m_mapParameterSet.contains(mhg_map))
		{
			VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
			if(!position.isValid())
			{
				CGoGNerr << "position attribute is not valid" << CGoGNendl;
				return;
			}

			MapParameters& mapParam = m_mapParameterSet[mhg_map];

			for(QHash<QString, Camera*>::iterator it = mapParam.depthCameraSet.begin(); it != mapParam.depthCameraSet.end(); ++it)
			{
				Camera* camera = it.value();
				QString cameraName(camera->getName());

				QString generatedName(mapName);
				generatedName.append("-");
				generatedName.append(cameraName);

				float near = camera->zNear();
				float far = camera->zFar();

				m_depthShader->setZMin(near/far);
				m_depthShader->setZMax(1.);

				m_depthShader->setAttributePosition(m_mapParameterSet[mhg_map].positionVBO);

				int width = m_depthFBO->getWidth(), height = m_depthFBO->getHeight();

				QImage image(width, height, QImage::Format_RGB32);

				m_schnapps->getSelectedView()->setCurrentCamera(camera);

				m_depthFBO->bind();
				glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );	//To get a clean texture (black background)
				mhg_map->draw(m_depthShader, CGoGN::Algo::Render::GL2::TRIANGLES);	//Render the map into the FrameBufferObject

				//Read pixels of the generated texture and store them in an image
				glReadPixels(0, 0, m_depthFBO->getWidth(), m_depthFBO->getHeight(), GL_RGBA, GL_UNSIGNED_BYTE, image.bits());
				m_depthFBO->unbind();

				image = image.mirrored();

				mapParam.depthImageSet[generatedName] = image;

				QString filename(directory);
				filename.append("/");
				filename.append(mapName);
				filename.append("/");

				mkdir(filename.toStdString().c_str(), 0777);

				filename.append(mapName);
				filename.append("-");
				filename.append(cameraName);
				filename.append(".png");

				if(!image.save(filename))
				{
					CGoGNerr << "Image '" << filename.toStdString() << "' has not been saved" << CGoGNendl;
				}

				MapHandlerGen* mhg_generated = m_schnapps->addMap(generatedName, 2);
				mapParam.projectedMapSet[generatedName] = mhg_generated;

				MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);
				PFP2::MAP* generated_map = mh_generated->getMap();

				VertexAttribute<PFP2::VEC3, PFP2::MAP> planeCoordinatesGenerated = mh_generated->getAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
				if(!planeCoordinatesGenerated.isValid())
				{
					planeCoordinatesGenerated = mh_generated->addAttribute<PFP2::VEC3, VERTEX>("PlaneCoordinates");
				}

				VertexAttribute<ImageCoordinates, PFP2::MAP> imageCoordinatesGenerated = mh_generated->getAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
				if(!imageCoordinatesGenerated.isValid())
				{
					imageCoordinatesGenerated = mh_generated->addAttribute<ImageCoordinates, VERTEX>("ImageCoordinates");
				}

				Algo::Surface::Tilings::Square::Grid<PFP2> grid(*generated_map, width-1, height-1);
				grid.embedIntoGrid(planeCoordinatesGenerated, width-1, height-1);

				std::vector<Dart> vDarts = grid.getVertexDarts();

				for(int i = 0; i < width; ++i)
				{
					for(int j = 0; j < height; ++j)
					{
						//Set planeCoordinates in [-1;1]
						planeCoordinatesGenerated[vDarts[j*width+i]][0] /= (width-1)/2.f;
						planeCoordinatesGenerated[vDarts[j*width+i]][1] /= (height-1)/2.f;

						imageCoordinatesGenerated[vDarts[j*width+i]].setCoordinates(i, height-j-1);
					}
				}

				mh_generated->notifyConnectivityModification();
				mh_generated->updateBB(planeCoordinatesGenerated);
				mh_generated->notifyAttributeModification(planeCoordinatesGenerated);
				mh_generated->notifyAttributeModification(imageCoordinatesGenerated);
				project2DImageTo3DSpace(mapName, generatedName);
			}
		}
	}
}

void Surface_DepthMapRendering_Plugin::project2DImageTo3DSpace(const QString& mapOrigin, const QString& mapGenerated)
{
	MapHandlerGen* mhg_origin = m_schnapps->getMap(mapOrigin);
	MapHandlerGen* mhg_generated = m_schnapps->getMap(mapGenerated);

	MapHandler<PFP2>* mh_origin = static_cast<MapHandler<PFP2>*>(mhg_origin);
	MapHandler<PFP2>* mh_generated = static_cast<MapHandler<PFP2>*>(mhg_generated);

	if(mh_origin && mh_generated && m_mapParameterSet.contains(mhg_origin))
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
		QImage image = mapParam.depthImageSet[mapGenerated];

		qglviewer::Vec camera_p = camera->position();
		PFP2::VEC3 camera_position = PFP2::VEC3(camera_p.x, camera_p.y, camera_p.z);

		PFP2::VEC3 plane_center_position = camera_position;
		plane_center_position[2] -= camera->zNear();

//		float d_C_Xc = (plane_center_position-camera_position).norm2();

		GLdouble p_matrix[16], mv_matrix[16];
		camera->getProjectionMatrix(p_matrix);
		camera->getModelViewMatrix(mv_matrix);

		PFP2::MATRIX44 projection_matrix, projection_matrix_inv, model_view_matrix, model_view_matrix_inv;

		for(int i = 0; i < 4; ++i)
		{
			for(int j = 0; j < 4; ++j)
			{
				projection_matrix(i,j) = p_matrix[i+4*j];
				model_view_matrix(i,j) = mv_matrix[i+4*j];
			}
		}

		float near = camera->zNear();
		float far = camera->zFar();

		projection_matrix.invert(projection_matrix_inv);
		model_view_matrix.invert(model_view_matrix_inv);

		TraversorV<PFP2::MAP> trav_vert_map(*generated_map);
		for(Dart d = trav_vert_map.begin(); d != trav_vert_map.end(); d = trav_vert_map.next())
		{
			float color = 1-qRed(image.pixel(imageCoordinatesGenerated[d].getXCoordinate(), imageCoordinatesGenerated[d].getYCoordinate()))/255.f;

			PFP2::VEC4 pos = PFP2::VEC4(planeCoordinatesGenerated[d][0], planeCoordinatesGenerated[d][1], color, 1.f);

			pos = model_view_matrix_inv*projection_matrix_inv*pos;

			positionGenerated[d] = PFP2::VEC3(pos[0]/pos[3], pos[1]/pos[3], pos[2]/pos[3]);
		}

		mh_generated->notifyAttributeModification(positionGenerated);
		mh_generated->updateBB(positionGenerated);
		m_schnapps->getSelectedView()->updateGL();
	}
}

#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_Plugin, Surface_DepthMapRendering_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_PluginD, Surface_DepthMapRendering_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
