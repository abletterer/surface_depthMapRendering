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

	m_depthFBO = new CGoGN::Utils::FBO(1024, 1024);
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

	if(mh_map)
	{
		PFP2::MAP* map = mh_map->getMap();

		QString base_name("DepthCamera-");

		m_cameraSet.reserve(number);

		for(int i = 0; i < number; ++i)
		{
			QString camera_name(base_name);
			camera_name.append(QString::number(i));
			Camera* camera = m_schnapps->addCamera(camera_name);

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
			m_cameraSet.push_back(camera);
		}
	}
}

void Surface_DepthMapRendering_Plugin::render(const QString& mapName, const QString& directory)
{
	if(m_depthShader)
	{
		MapHandlerGen* mhg_map = m_schnapps->getMap(mapName);
		MapHandler<PFP2>* mh_map = static_cast<MapHandler<PFP2>*>(mhg_map);

		if(mh_map)
		{
			VertexAttribute<PFP2::VEC3, PFP2::MAP> position = mh_map->getAttribute<PFP2::VEC3, VERTEX>("position");
			if(!position.isValid())
			{
				CGoGNerr << "position attribute is not valid" << CGoGNendl;
				return;
			}

			for(std::vector<Camera*>::const_iterator cam = m_cameraSet.begin(); cam != m_cameraSet.end(); ++cam)
			{
				float near = (*cam)->zNear();
				float far = (*cam)->zFar();

				m_depthShader->setZMin(near/far);
				m_depthShader->setZMax(1.);

				m_depthShader->setAttributePosition(m_mapParameterSet[mhg_map].positionVBO);

				int width = m_depthFBO->getWidth(), height = m_depthFBO->getHeight();

				QImage image(width, height, QImage::Format_RGB32);

				m_schnapps->getSelectedView()->setCurrentCamera(*cam);

				m_depthFBO->bind();
				glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );	//To get a clean texture (black background)
				mhg_map->draw(m_depthShader, CGoGN::Algo::Render::GL2::TRIANGLES);	//Render the map into the FrameBufferObject

				//Read pixels of the generated texture and store them in a vector
				glReadPixels(0, 0, m_depthFBO->getWidth(), m_depthFBO->getHeight(), GL_RGBA, GL_UNSIGNED_BYTE, image.bits());
				m_depthFBO->unbind();

				image = image.mirrored();

				QString filename(directory);
				filename.append("/");
				filename.append(mapName);
				filename.append("/");

				mkdir(filename.toStdString().c_str(), 0777);

				filename.append(mapName);
				filename.append("-");
				filename.append((*cam)->getName());
				filename.append(".png");

				if(!image.save(filename))
				{
					CGoGNerr << "Image '" << filename.toStdString() << "' has not been saved" << CGoGNendl;
				}
			}
		}
	}
}

#ifndef DEBUG
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_Plugin, Surface_DepthMapRendering_Plugin)
#else
Q_EXPORT_PLUGIN2(Surface_DepthMapRendering_PluginD, Surface_DepthMapRendering_Plugin)
#endif

} // namespace SCHNApps

} // namespace CGoGN
