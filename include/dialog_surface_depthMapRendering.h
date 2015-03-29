#ifndef _DIALOG_SURFACE_DEPTHMAPRENDERING_H_
#define _DIALOG_SURFACE_DEPTHMAPRENDERING_H_

#include "ui_dialog_surface_depthMapRendering.h"

namespace CGoGN
{

namespace SCHNApps
{

class SCHNApps;
class MapHandlerGen;
class Surface_DepthMapRendering_Plugin;

class Dialog_Surface_DepthMapRendering: public QDialog, public Ui::Dialog_Surface_DepthMapRendering
{
	Q_OBJECT


	friend class Surface_DepthMapRendering_Plugin;

public:
	Dialog_Surface_DepthMapRendering(SCHNApps* s);

private:

	void updateMapParameters();

private:
	SCHNApps* m_schnapps;
	MapHandlerGen* m_selectedMap;

public slots:
	void selectedMapChanged();
	void addMapToList(MapHandlerGen* map);
	void removeMapFromList(MapHandlerGen* map);
};

} // namespace SCHNApps

} // namespace CGoGN

#endif
