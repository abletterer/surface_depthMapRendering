#ifndef _IMAGECOORDINATES_H_
#define _IMAGECOORDINATES_H_

namespace CGoGN
{
namespace SCHNApps
{

/*
* Classe définissant l'attribut de sommet Coordinates, définissant les coordonnées images associées à un point
*/
class ImageCoordinates
{
public:
	ImageCoordinates()
		: m_x(0), m_y(0)
	{}

	~ImageCoordinates()
	{
	}

	void setCoordinates(int x, int y)
	{
		m_x = x;
		m_y = y;
	}

	int getXCoordinate() { return m_x; }
	void setXCoordinate(int x)
	{
		m_x = x;
	}

	int getYCoordinate() { return m_y; }
	void setYCoordinate(int y)
	{
		m_y = y;
	}

	friend std::ostream& operator<< (std::ostream& stream, const ImageCoordinates &coordinates)
	{
		return stream;
	}

	static std::string CGoGNnameOfType()
	{
		return "ImageCoordinates" ;
	}

private:
	int m_x, m_y;
};

} //namespace SCHNApps
} //namespace CGoGN

#endif
