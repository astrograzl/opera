#ifndef OPERAGEOMETRICSHAPES_H
#define OPERAGEOMETRICSHAPES_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaGeometricShapes
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Apr/2012
 
 Copyright (C) 2011  Opera Pipeline team, Canada France Hawaii Telescope
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see:
 http://software.cfht.hawaii.edu/licenses
 -or-
 http://www.gnu.org/licenses/gpl-3.0.html
 ********************************************************************/

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$


#include <ostream>
#include <math.h>

#include "libraries/operaLibCommon.h"			// for MAXNPOLYGONSIDES
#include "libraries/operaStats.h"   // for operaArrayIndexSort

class operaPoint;
class Rectangle;
class Circle;
class Polygon;
class Line;

/*! 
 * \sa class operaGeometricShapes
 * \brief geometric shapes containers.
 * \details The operaGeometricShapes defines geometric shapes that can be used  
 * \details by operaExtractionAperture in order to obtain information from 
 * \details specific regions of FITS images.
 * \return none
 * \file operaGeometricShapes.h
 * \ingroup libraries
 */

/*
 * class operaPoint
 */
class operaPoint {
private:
    float Xcoord, Ycoord;
    
public:
    operaPoint(void);
    
    operaPoint(float xp, float yp);    
    
    operaPoint(Line &inputLine1, Line &inputLine2);
    
	operaPoint(operaPoint &point);
    
    ~operaPoint(void);
    
	void setPoint(operaPoint &point);
    
    void setPoint(float xp, float yp);

    float getXcoord(void);
    
    float getYcoord(void);
    
};

/*
 * Simple boxes for things like guide windows in WIRCam
 */
class Box {
	
private:
	unsigned x1, x2, y1, y2;
	unsigned xDim, yDim;
	
public:
	
	Box();
	
	Box(unsigned X1, unsigned X2, unsigned Y1, unsigned Y2);
	
	Box(unsigned X1, unsigned X2, unsigned Y1, unsigned Y2, unsigned XDim, unsigned YDim);
	
	unsigned getX1(void);
	
	unsigned getX2(void);
	
	unsigned getY1(void);
	
	unsigned getY2(void);
	
	unsigned getDX(void);	
	
	unsigned getDY(void);
	
	unsigned getXDim(void);
	
	unsigned getYDim(void);
	
	void setX1(unsigned X1);
	
	void setX2(unsigned X2);
	
	void setY1(unsigned Y1);
	
	void setY2(unsigned Y2);
	
	void setXDim(unsigned XDim);
	
	void setYDim(unsigned YDim);
	
	unsigned getSize(void);
	
};

/*
 * class Rectangle
 */
class Rectangle {
private:
    float width;
    float height;
    float angle; // 0 <= angle < 180 degrees
    operaPoint center;
    operaPoint corners[FOURSIDES];
    
public:
    Rectangle(void);
    
    Rectangle(float Width, float Height, float Angle);
	
    Rectangle(float Width, float Height, float Angle, operaPoint &Center);   
    
    Rectangle(float Width, float Height, float Angle, operaPoint *Center);   
    
	Rectangle(Rectangle &Rectangle);
    
    ~Rectangle(void);
    
	void setRectangle(float Width, float Height, float Angle, operaPoint &Center);
    
	void setRectangle(Rectangle &Rectangle);
    
	operaPoint* getCorner(unsigned index);
    
    float getWidth(void);
    
    float getHeight(void);
    
    float getAngle(void);   
	
    operaPoint *getCenter(void);
    
    void printCorners(void); 
	
    bool pointInRectangle(operaPoint &TestPoint);   
};

/*
 * class Circle
 */
class Circle {
private:
    float radius;
    operaPoint center;   
    
public:
    Circle(void);     
	
    Circle(float Radius);
    
    Circle(float Radius, operaPoint &Center);
    
	Circle(Circle &ACircle);
    
    ~Circle(void);
    
	void setCircle(Circle &ACircle);
    
    float getRadius(void);
    
    operaPoint *getCenter(void);
    
    bool pointInCircle(operaPoint &TestPoint);
};

/*
 * class Polygon
 */
class Polygon {
private:
    unsigned nSides;
    operaPoint vertices[MAXNPOLYGONSIDES];
    
public:
    Polygon(void);
    
    Polygon(unsigned NSides, operaPoint *Vertices[]);
    
	Polygon(unsigned NSides, operaPoint Vertices[]);
 	
	Polygon(Polygon &Poly);
	
    ~Polygon(void);
	
	void setPolygon(Polygon &Poly);
	
    void simplePolygonization(void);
    
    void setNSides(unsigned NSides);
    
    unsigned getNSides(void);
    
    operaPoint* getVertex(unsigned index);
    
    bool pointInPolygon(operaPoint &testPoint);
    
    void printVertexCoordinates(void);    
    
};

/*
 * class Line
 */
typedef enum LinePosition {
	line_duplicate, line_top, line_left, line_right, line_bottom, line_perpendicular
} LinePosition_t;

class Line {
private:
    float slope;
    float intercept;    
    float width;
    float length;
    operaPoint midPoint;
    
public:
    Line(void);
    
    Line(float Slope);
    
    Line(float Slope, operaPoint &SamplePoint);
    
    Line(float Slope, float Intercept);
    
    Line(float Slope, float Width, float Length);
    
    Line(float Slope, float Width, float Length, operaPoint &MidPoint);
    
    Line(float Slope, float Width, float Length, operaPoint *MidPoint);
    
    Line(float Slope, float Intercept, float Width, float Length);
	
	Line(Line &inputLine, LinePosition_t LinePosition);
	
	Line(Line &ALine);
    
    ~Line(void);
    
	void setLine(Line &ALine);
    
    float getSlope(void);
    
    float getWidth(void);
    
    float getLength(void);   
    
    operaPoint *getMidPoint(void);
    
    void setMidPoint(operaPoint &MidPoint);    
	
    void setMidPoint(operaPoint *MidPoint);    
	
    float getIntercept(void);
    
    void printLineEquation(void); 
    
    float getYcoord(float x); 
    
    float getXcoord(float y); 
    
    bool pointOnLine(operaPoint &TestPoint);
    
    bool pointInLineWidth(operaPoint &TestPoint, Line &inputLine);    
    
    Line *newPerpendicularLine(Line &inputLine);  
    
	operaPoint *newIntersectionPoint(Line &inputLine);
    
    Line *newTopLine(Line &inputLine);
    
    Line *newBottomLine(Line &inputLine);
    
    Line *newLeftLine(Line &inputLine);
    
    Line *newRightLine(Line &inputLine); 
    
    float getLineYWidth(void);     
};

/*
 class Donut {
 private:
 float innerRadius;
 float outerRadius;
 operaPoint center;
 public:
 
 }
 */

#endif

