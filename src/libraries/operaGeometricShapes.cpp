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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaGeometricShapes.h"
#include "libraries/operaException.h"

#include "libraries/operaLibCommon.h"
#include "libraries/operaMatrix.h"

/*!
 * operaGeometricShapes
 * \author Doug Teeple / Eder Martioli
 * \brief GeometricShapes  
 * \details The operaGeometricShapes defines geometric shapes that can be used  
 * \details by operaExtractionAperture in order to obtain information from 
 * \details specific regions of FITS images.
 * \file operaGeometricShapes.cpp
 * \ingroup libraries
 */

using namespace std;


/*
 * Constructors
 */


/*
 * operaPoint class 
 */

operaPoint::operaPoint(void) :
Xcoord(0),
Ycoord(0)
{
}

operaPoint::operaPoint(float xp, float yp) {
    Xcoord = xp;
    Ycoord = yp;
}

operaPoint::operaPoint(const Line &inputLine1, const Line &inputLine2) {
    if(inputLine1.getSlope() == inputLine2.getSlope()) {
        throw operaException("operaPoint: error: no intersection point; input lines are parallel. ", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);                
    }
    
    float xintersec = (inputLine1.getIntercept() - inputLine2.getIntercept())/(inputLine2.getSlope()-inputLine1.getSlope());
    float yintersec = inputLine1.getYcoord(xintersec);
    
    Xcoord = xintersec;
    Ycoord = yintersec;    
}

operaPoint::operaPoint(const operaPoint &APoint) {
    Xcoord = APoint.Xcoord;
    Ycoord = APoint.Ycoord;
}

void operaPoint::setPoint(const operaPoint &APoint) {
    Xcoord = APoint.Xcoord;
    Ycoord = APoint.Ycoord;
}

/*
 * Simple boxes for things like guide windows in WIRCam
 */
Box::Box() :
x1(0),
x2(0),
y1(0),
y2(0),
xDim(0),
yDim(0)
{
}

Box::Box(unsigned X1, unsigned X2, unsigned Y1, unsigned Y2) :
xDim(0),
yDim(0)
{
	x1 = X1;
	x2 = X2;
	y1 = Y1;
	y2 = Y2;
}

Box::Box(unsigned X1, unsigned X2, unsigned Y1, unsigned Y2, unsigned XDim, unsigned YDim)
{
	x1 = X1;
	x2 = X2;
	y1 = Y1;
	y2 = Y2;
	xDim = XDim;
	yDim = YDim;
}

unsigned Box::getX1(void) const {
	return x1;
}

unsigned Box::getX2(void) const {
	return x2;
}

unsigned Box::getY1(void) const {
	return y1;
}

unsigned Box::getY2(void) const {
	return y2;
}

unsigned Box::getDX(void) const {
	return x2-x1;
}

unsigned Box::getDY(void) const {
	return y2-y1;
}

void Box::setX1(unsigned X1) {
	x1 = X1;
}

void Box::setX2(unsigned X2) {
	x2 = X2;
}

void Box::setY1(unsigned Y1) {
	y1 = Y1;
}

void Box::setY2(unsigned Y2) {
	y2 = Y2;
}

unsigned Box::getSize(void) const {
	return ((y2-y1)*(x2-x1));
}

/*
 * Rectangle class 
 */

Rectangle::Rectangle(void) :
width(0),
height(0),
angle(0),
center(0,0) 
{
}

Rectangle::Rectangle(float Width, float Height, float Angle) :
center(0,0) 
{   
    if(Width <= 0 || Height <= 0){
        throw operaException("Rectangle: error: sides must be greater than zero",operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }   
    width = Width;
    height = Height;
    angle = Angle;
    
    operaPoint corners_tmp[FOURSIDES];
    
    corners_tmp[0].setPoint(-width/2,-height/2);
    corners_tmp[1].setPoint(+width/2,-height/2);    
    corners_tmp[2].setPoint(+width/2,+height/2);    
    corners_tmp[3].setPoint(-width/2,+height/2);    
    
    CMatrix rotationMatrix = newCMatrix(2,2);
    RotationMatrix2D(angle, rotationMatrix);
    CMatrix pointMatrix = newCMatrix(1,2);
    CMatrix rotatedPointMatrix = newCMatrix(1,2);
    
    for(unsigned i=0;i<FOURSIDES;i++) {
        pointMatrix[0][0] = corners_tmp[i].getXcoord();
        pointMatrix[1][0] = corners_tmp[i].getYcoord();        
        MatrixMultiplication(rotationMatrix,pointMatrix, rotatedPointMatrix);    
        corners[i].setPoint(rotatedPointMatrix[0][0],rotatedPointMatrix[1][0]);
    }
	deleteCMatrix(rotatedPointMatrix);
	deleteCMatrix(rotationMatrix);
    deleteCMatrix(pointMatrix);
}

Rectangle::Rectangle(float Width, float Height, float Angle, const operaPoint &Center) {   
    if(Width <= 0 || Height <= 0){
        throw operaException("Rectangle: error: sides must be greater than zero",operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }   
    width = Width;
    height = Height;
    angle = Angle;
    center.setPoint(Center);
    
    operaPoint corners_tmp[FOURSIDES];
    
    corners_tmp[0].setPoint(-width/2,-height/2);
    corners_tmp[1].setPoint(+width/2,-height/2);    
    corners_tmp[2].setPoint(+width/2,+height/2);    
    corners_tmp[3].setPoint(-width/2,+height/2);    
	
    CMatrix rotationMatrix = newCMatrix(2,2);
    RotationMatrix2D(angle, rotationMatrix);
    CMatrix pointMatrix = newCMatrix(1,2);
    CMatrix rotatedPointMatrix = newCMatrix(1,2);
    
    for(unsigned i=0;i<FOURSIDES;i++) {
        pointMatrix[0][0] = corners_tmp[i].getXcoord();
        pointMatrix[1][0] = corners_tmp[i].getYcoord();
        MatrixMultiplication(rotationMatrix,pointMatrix, rotatedPointMatrix);
        corners[i].setPoint(rotatedPointMatrix[0][0] + center.getXcoord(),rotatedPointMatrix[1][0] + center.getYcoord());
    }
	deleteCMatrix(rotatedPointMatrix);       
    deleteCMatrix(rotationMatrix);
    deleteCMatrix(pointMatrix);
}

Rectangle::Rectangle(float Width, float Height, float Angle, const operaPoint *Center) {   
    if(Width <= 0 || Height <= 0){
        throw operaException("Rectangle: error: sides must be greater than zero",operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }   
    width = Width;
    height = Height;
    angle = Angle;
    center.setPoint(*Center);
    
    operaPoint corners_tmp[FOURSIDES];
    
    corners_tmp[0].setPoint(-width/2,-height/2);
    corners_tmp[1].setPoint(+width/2,-height/2);    
    corners_tmp[2].setPoint(+width/2,+height/2);    
    corners_tmp[3].setPoint(-width/2,+height/2);    
    
    CMatrix rotationMatrix = newCMatrix(2,2);
    rotationMatrix = RotationMatrix2D(angle, rotationMatrix);
    CMatrix pointMatrix = newCMatrix(1,2);
    CMatrix rotatedPointMatrix = newCMatrix(1,2);
    
    for(unsigned i=0;i<FOURSIDES;i++) {
        pointMatrix[0][0] = corners_tmp[i].getXcoord();
        pointMatrix[1][0] = corners_tmp[i].getYcoord();
        MatrixMultiplication(rotationMatrix,pointMatrix,rotatedPointMatrix);
        corners[i].setPoint(rotatedPointMatrix[0][0] + center.getXcoord(),rotatedPointMatrix[1][0] + center.getYcoord());
    }
	deleteCMatrix(rotatedPointMatrix);       
    deleteCMatrix(rotationMatrix);
    deleteCMatrix(pointMatrix);
}

Rectangle::Rectangle(const Rectangle &aRectangle) {   
    width = aRectangle.width;
    height = aRectangle.height;
    angle = aRectangle.angle;
    center.setPoint(aRectangle.center);
    
    for(unsigned i=0;i<FOURSIDES;i++) {
        corners[i].setPoint(aRectangle.corners[i]);
    }
}

void Rectangle::setRectangle(const Rectangle &aRectangle) {
    width = aRectangle.width;
    height = aRectangle.height;
    angle = aRectangle.angle;
    center.setPoint(aRectangle.center);
    
    for(unsigned i=0;i<FOURSIDES;i++) {
        corners[i].setPoint(aRectangle.corners[i]);
    }
}

void Rectangle::setRectangle(float Width, float Height, float Angle, const operaPoint &Center) {
    width = Width;
    height = Height;
    angle = Angle;
    center.setPoint(Center);
    
    operaPoint corners_tmp[FOURSIDES];
    
    corners_tmp[0].setPoint(-width/2,-height/2);
    corners_tmp[1].setPoint(+width/2,-height/2);    
    corners_tmp[2].setPoint(+width/2,+height/2);    
    corners_tmp[3].setPoint(-width/2,+height/2);    
    
    CMatrix rotationMatrix = newCMatrix(2,2);
    rotationMatrix = RotationMatrix2D(angle, rotationMatrix);
    CMatrix pointMatrix = newCMatrix(1,2);
    CMatrix rotatedPointMatrix = newCMatrix(1,2);
    
    for(unsigned i=0;i<FOURSIDES;i++) {
        pointMatrix[0][0] = corners_tmp[i].getXcoord();
        pointMatrix[1][0] = corners_tmp[i].getYcoord();
        rotatedPointMatrix =  MatrixMultiplication(rotationMatrix,pointMatrix, rotatedPointMatrix);
        corners[i].setPoint(rotatedPointMatrix[0][0] + center.getXcoord(),rotatedPointMatrix[1][0] + center.getYcoord());
    }
	deleteCMatrix(rotatedPointMatrix);    
    deleteCMatrix(rotationMatrix);
    deleteCMatrix(pointMatrix);
}

/*
 * Circle class 
 */

Circle::Circle(void)  :
radius(0),
center(0,0)
{
}

Circle::Circle(float Radius)  :
center(0,0)
{
    if(Radius <= 0){
        throw operaException("Circle: error: radius must be greater than zero",operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }          
    radius = Radius;
}

Circle::Circle(float Radius, const operaPoint &Center) {
    if(Radius <= 0){
        throw operaException("Circle: error: radius must be greater than zero",operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }          
    radius = Radius;
    center = Center;   
}

Circle::Circle(const Circle &ACircle) {
    radius = ACircle.radius;
    center.setPoint(ACircle.center);   
}

void Circle::setCircle(const Circle &ACircle) {
    radius = ACircle.radius;
    center.setPoint(ACircle.center);   
}

/*
 * Polygon class 
 */

Polygon::Polygon(void) :
nSides(0)
{ 
}

Polygon::Polygon(unsigned NSides, const operaPoint *Vertices[]) {
    if(NSides < 3){
        throw operaException("Polygon: error: number of sides must be 3 or greater", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    nSides = NSides;
    
    for(unsigned i=0; i<nSides; i++) {
        vertices[i].setPoint((*Vertices)[i]);
    }
    
    simplePolygonization(); // organize vertices to make polygon simple, i.e. with no intersecting lines     
}

Polygon::Polygon(unsigned NSides, const operaPoint Vertices[]) {
    if(NSides < 3){
        throw operaException("Polygon: error: number of sides must be 3 or greater", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    nSides = NSides;
    
    for(unsigned i=0; i<nSides; i++) {
        vertices[i].setPoint(Vertices[i]);
    }
    
    simplePolygonization(); // organize vertices to make polygon simple, i.e. with no intersecting lines     
}

Polygon::Polygon(const Polygon &Poly) {
    
    nSides = Poly.nSides;
    
    for(unsigned i=0; i<nSides; i++) {
        vertices[i].setPoint(Poly.vertices[i]);
    }
    
    simplePolygonization(); // organize vertices to make polygon simple, i.e. with no intersecting lines     
}

void Polygon::setPolygon(const Polygon &Poly) {
    
    nSides = Poly.nSides;
    
    for(unsigned i=0; i<nSides; i++) {
        vertices[i].setPoint(Poly.vertices[i]);
    }
    
    simplePolygonization(); // organize vertices to make polygon simple, i.e. with no intersecting lines     
}

/*
 * Line class 
 */

Line::Line(void) :
slope(0),
intercept(0),
width(0),
length(BIG),
midPoint(0,0)
{ 
}

Line::Line(float Slope) :
intercept(0),
width(0),
length(BIG),
midPoint(0,0)
{
    slope = Slope;   
}

Line::Line(float Slope, float Intercept) :
width(0),
length(BIG)
{
    slope = Slope;   
    intercept = Intercept;
	midPoint.setPoint(0, intercept);
}

Line::Line(float Slope, const operaPoint &SamplePoint) :
width(0),
length(BIG)
{
    slope = Slope;
	midPoint.setPoint(SamplePoint.getXcoord(), SamplePoint.getYcoord());
    intercept = midPoint.getYcoord() - slope*midPoint.getXcoord();
}

Line::Line(float Slope, float Width, float Length) :
intercept(0),
midPoint(0,0)
{
    if(Width < 0 || Length <= 0){
        throw operaException("Line: error: line dimensions must be greater than zero", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }          
    slope = Slope;
    width = Width;
    length = Length;
}

Line::Line(float Slope, float Width, float Length, const operaPoint &MidPoint) {
    if(Width < 0 || Length <= 0){
        throw operaException("Line: error: line dimensions must be greater than zero", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }          
    slope = Slope;
    width = Width;
    length = Length;
	setMidPoint(MidPoint);
    intercept = midPoint.getYcoord() - slope*midPoint.getXcoord();
}

Line::Line(float Slope, float Width, float Length, const operaPoint *MidPoint) {
    if(Width < 0 || Length <= 0){
        throw operaException("Line: error: line dimensions must be greater than zero", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }          
    slope = Slope;
    width = Width;
    length = Length;
	setMidPoint(MidPoint);
    intercept = midPoint.getYcoord() - slope*midPoint.getXcoord();
}

Line::Line(float Slope, float Intercept, float Width, float Length) {
    if(Width < 0 || Length <= 0){
        throw operaException("Line: error: line dimensions must be greater than zero", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }          
    slope = Slope;
    width = Width;
    length = Length;
    intercept = Intercept;
	midPoint.setPoint(0, intercept);
}

Line::Line(const Line &inputLine) {
	slope = inputLine.slope;
	width = inputLine.width;
	length = inputLine.length;
	intercept = inputLine.intercept;
	setMidPoint(inputLine.midPoint); 
}

/*
 * Create a line of line position type LinePosition_t
 */
Line::Line(const Line &inputLine, LinePosition_t LinePosition)
{
	switch (LinePosition) {
		case line_duplicate: {
			slope = inputLine.slope;
			width = inputLine.width;
			length = inputLine.length;
			intercept = inputLine.intercept;
			setMidPoint(inputLine.midPoint); 
		}
			break;
		case line_top: {
			setLine(inputLine);
			float ShiftInYDirection = width/(2.0*cos(atan(slope)));
			Line perpendicularline(inputLine, line_perpendicular);
			Line topline(slope,(intercept+ShiftInYDirection),width,length);
			operaPoint midPoint(topline, perpendicularline);
			slope = topline.slope;
			width = topline.width;
			length = topline.length;
			intercept = topline.intercept;
			setMidPoint(midPoint); 
		}
			break;
		case line_bottom: {
			setLine(inputLine);
			float ShiftInYDirection = width/(2.0*cos(atan(slope)));
			Line perpendicularline(inputLine, line_perpendicular);
			Line bottomline(slope,(intercept-ShiftInYDirection),width,length);
			operaPoint midPoint(bottomline, perpendicularline);
			slope = bottomline.slope;
			width = bottomline.width;
			length = bottomline.length;
			intercept = bottomline.intercept;
			setMidPoint(midPoint); 
		}
			break;
		case line_left: {
			setLine(inputLine);
			float midPointXshift = (length/2.0)*cos(atan(slope));
			Line perpendicularline(inputLine, line_perpendicular);
			operaPoint leftMidPoint(midPoint.getXcoord() - midPointXshift, getYcoord(midPoint.getXcoord() - midPointXshift));
			Line leftline(perpendicularline.slope,length,width,leftMidPoint);
			slope = leftline.slope;
			width = leftline.width;
			length = leftline.length;
			intercept = leftline.intercept;
			setMidPoint(leftMidPoint); 
		}
			break;
		case line_right: {
			setLine(inputLine);
			float midPointXshift = (length/2.0)*cos(atan(slope));
			Line perpendicularline(inputLine, line_perpendicular);
			operaPoint rightMidPoint(midPoint.getXcoord() + midPointXshift, getYcoord(midPoint.getXcoord() + midPointXshift));
			Line rightline(perpendicularline.slope,length,width,rightMidPoint);
			slope = perpendicularline.slope;
			width = rightline.width;
			length = rightline.length;
			intercept = rightline.intercept;
			setMidPoint(rightMidPoint); 
		}
			break;
		case line_perpendicular: {
			slope = -(1.0/inputLine.slope);
			width = inputLine.length;
			length = inputLine.width;
			setMidPoint(inputLine.getMidPoint());
			intercept = midPoint.getYcoord() - slope*midPoint.getXcoord();
		}
			break;
		default:
			break;
	}
}

void Line::setLine(const Line &ALine) {
    slope = ALine.slope;
    width = ALine.width;
    length = ALine.length;
    intercept = ALine.intercept;
	midPoint.setPoint(ALine.midPoint.getXcoord(), ALine.midPoint.getYcoord());
}

/*
 * Destructors
 */

operaPoint::~operaPoint(void) {
    
}

Rectangle::~Rectangle(void) {
}

Circle::~Circle(void) {
}

Polygon::~Polygon(void) {
    nSides = 0;
}

Line::~Line(void) {
}


/*
 * Setter/Getters
 */

float operaPoint::getXcoord(void) const {
    return Xcoord;
}

float operaPoint::getYcoord(void) const {
    return Ycoord;
}

void operaPoint::setPoint(float xp, float yp) {
    Xcoord = xp;
    Ycoord = yp;
}

float Circle::getRadius(void) const {
    return radius;
}

operaPoint* Circle::getCenter(void) {
    return &center;
}

const operaPoint* Circle::getCenter(void) const {
    return &center;
}

operaPoint* Rectangle::getCorner(unsigned index) {
    return &corners[index];
}

const operaPoint* Rectangle::getCorner(unsigned index) const {
    return &corners[index];
}

float Rectangle::getWidth(void) const {
	return width;   
}

float Rectangle::getHeight(void) const {
    return height;    
}

float Rectangle::getAngle(void) const {
    return angle;    
}

operaPoint* Rectangle::getCenter(void) {
    return &center;
}

const operaPoint* Rectangle::getCenter(void) const {
    return &center;
}

float Line::getSlope(void) const {
	return slope;     
}

float Line::getWidth(void) const {
	return width;     
}

float Line::getLength(void) const {
	return length;     
}

operaPoint *Line::getMidPoint(void) {
	return &midPoint;     
}

const operaPoint *Line::getMidPoint(void) const {
	return &midPoint;     
}

void Polygon::setNSides(unsigned NSides){
    nSides = NSides;   
}

unsigned Polygon::getNSides(void) const {
    return nSides;     
}

operaPoint* Polygon::getVertex(unsigned index){
    return &vertices[index];     
}

const operaPoint* Polygon::getVertex(unsigned index) const {
    return &vertices[index];     
}

/*
 * Methods
 */

void Polygon::simplePolygonization(void) {
    unsigned anchorPoint=0;
    float minYCoord = 1e30;
    
    // First determine the anchor point as the point with minimum Y coordinate value. 
    for(unsigned i=0;i<nSides;i++){
        if(vertices[i].getYcoord() < minYCoord) {
            minYCoord = vertices[i].getYcoord();
            anchorPoint = i;            
        }
    } 
	
    float polarAngle[MAXNPOLYGONSIDES];
    
    unsigned np=0;
    for(unsigned i=0;i<nSides;i++){
        if(i != anchorPoint) {
            float x = vertices[i].getXcoord() - vertices[anchorPoint].getXcoord();
            float y = vertices[i].getYcoord() - vertices[anchorPoint].getYcoord();
            polarAngle[np] = atan2(y,x);
			
#ifdef PRINT_DEBUG    
            cerr << "Polygon::simplePolygonization: anchorPoint = " << anchorPoint << endl;
            cerr << "Polygon::simplePolygonization: i = "<< i << " vertexXcoords,vertexYcoords = " << vertices[i]->getXcoord() << "," << vertices[i]->getYcoord() << endl;
            cerr << "Polygon::simplePolygonization: x,y = " << x << "," << y << endl;
            cerr << "Polygon::simplePolygonization: np = " << np << " polarAngle = " <<  polarAngle[np] << endl;
#endif            
            np++;
        }
    } 
    int sindex[MAXNPOLYGONSIDES];
    operaArrayIndexSort((int)np,polarAngle,sindex);
    
    operaPoint NewVertices[MAXNPOLYGONSIDES]; 
    
    NewVertices[0].setPoint(vertices[anchorPoint].getXcoord(), vertices[anchorPoint].getYcoord());
    
    np=0;
    for(unsigned i=0;i<nSides;i++){
        if(i != anchorPoint) {   
            NewVertices[np+1].setPoint(vertices[sindex[np]].getXcoord(), vertices[sindex[np]].getYcoord());            
            np++;
        }            
    }
    
    for(unsigned i=0;i<nSides;i++){
        vertices[i].setPoint(NewVertices[i]);  
    }
}


/* 
 * Polygon::pointInPolygon(operaPoint testPoint)
 * \author Doug Teeple / Eder Martioli
 * \brief This function will return TRUE if the test point is inside the polygon, or
 * \brief FALSE if it is not.  If the point is exactly on the edge of the polygon,
 * \brief then the function may return TRUE or FALSE.
 * \Notes Source: http://alienryderflex.com/polygon/
 * \ingroup libraries
 */

bool Polygon::pointInPolygon(const operaPoint &testPoint) const {
    
    // Comment: Eder, April 10 2012. 
    // This function still doesn't work properly
	// Eder addedparens for ambiguity - please check that this is what you want...
	// DT changed ^= (bitwise XOR) to |= (or equals)
    
    unsigned j=nSides-1 ;
    bool oddNodes = false;
    
    float testXcoord = testPoint.getXcoord();
    float testYcoord = testPoint.getYcoord();
    
    for (unsigned i=0; i<nSides; i++) {
        if (((vertices[i].getYcoord() < testYcoord && vertices[j].getYcoord()>=testYcoord)
             ||   (vertices[j].getYcoord() < testYcoord && vertices[i].getYcoord()>=testYcoord))
            &&  (vertices[i].getXcoord()<=testXcoord || vertices[j].getXcoord()<=testXcoord)) {
            
            oddNodes |= (vertices[i].getXcoord()+(testYcoord-vertices[i].getYcoord())/(vertices[j].getYcoord()-vertices[i].getYcoord())*(vertices[j].getXcoord()-vertices[i].getXcoord())) < testXcoord;
        }
        j=i; 
    }
    
    return oddNodes; 
}

void Polygon::printVertexCoordinates(void) {
    for (unsigned i=0; i<nSides; i++) {
        cout << "(" << vertices[i].getXcoord() << "," << vertices[i].getYcoord() << ")\t";
    }
    cout << endl;
}

bool Circle::pointInCircle(const operaPoint &TestPoint) const {
	
    bool PointIsInside = FALSE;
    
    float testXcoord = TestPoint.getXcoord();
    float testYcoord = TestPoint.getYcoord();
	
    float xsqr = (testXcoord - center.getXcoord())*(testXcoord - center.getXcoord());
    
    float yplus = center.getYcoord() + sqrt(radius*radius - xsqr);
    float yminus = center.getYcoord() - sqrt(radius*radius - xsqr);
    
    if(testYcoord > yminus && testYcoord < yplus) {
        PointIsInside = TRUE;
    }
    
    return PointIsInside;  
}

void Rectangle::printCorners(void) {
    for (unsigned i=0; i<FOURSIDES; i++) {
        cout << "(" << corners[i].getXcoord() << "," << corners[i].getYcoord() << ")\t";
    }
    cout << endl;    
}

bool Rectangle::pointInRectangle(const operaPoint &TestPoint) const {
	
    float lineSlope;
    float lineWidth;
    float lineLength;
	
    if(getAngle() == 90 || getAngle() == -90 || getAngle() == 270  || getAngle() == -270) {
        lineSlope = 0;
        lineWidth = getWidth();
        lineLength = getHeight();        
    } else {
        lineSlope = tan(getAngle()*M_PI/180.0);
        lineWidth = getHeight();
        lineLength = getWidth();
    }
	
    const operaPoint* lineMidPoint = getCenter();
    
    Line rectangleCenterLine(lineSlope,lineWidth,lineLength,lineMidPoint); 
    
    bool PointIsIN = rectangleCenterLine.pointInLineWidth(TestPoint, rectangleCenterLine);
    
    return PointIsIN;
}

void Line::printLineEquation(void) {
    cout << "y = " << slope << "*x + " << intercept << endl;
}

float Line::getYcoord(float x) const {
    return slope*x + intercept;
}

float Line::getXcoord(float y) const {
    if(slope == 0) {
        throw operaException("Line: error: x is undetermined because slope is zero", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);        
    }
	
    return (y - intercept)/slope;
}

/*
 * These need to be fixed. They generate new points each time -- that the caller has to clean up.
 */
Line *Line::newPerpendicularLine(const Line &inputLine) const {
    
    float perpendicularSlope = -(1./inputLine.getSlope());
    
    float perpendicularWidth = inputLine.getLength();
    
    float perpendicularLength = inputLine.getWidth();
    
    Line *perpendicularLine = new Line(perpendicularSlope,perpendicularWidth,perpendicularLength,inputLine.getMidPoint());
    
    return perpendicularLine;
}

Line *Line::newTopLine(const Line &inputLine) const {
    float ShiftInYDirection = inputLine.getLineYWidth()/2.0;
    
    Line *perpendicularline = newPerpendicularLine(inputLine);
    
    Line *topline = new Line(slope,(intercept+ShiftInYDirection),width,length); 
    
    operaPoint topMidPoint(*topline,*perpendicularline);
    
    topline->setMidPoint(topMidPoint); 
	
    delete perpendicularline;    
    
    return topline;
}

Line *Line::newBottomLine(const Line &inputLine) const {
    float ShiftInYDirection = inputLine.getLineYWidth()/2.0;
    
    Line *perpendicularline = newPerpendicularLine(inputLine);
    
    Line *bottomline = new Line(slope,(intercept-ShiftInYDirection),width,length);
    
    operaPoint bottomMidPoint(*bottomline,*perpendicularline);
    
    bottomline->setMidPoint(bottomMidPoint);
	
    delete perpendicularline;    
	
    return bottomline;
}

Line *Line::newLeftLine(const Line &inputLine) const {
    float midPointXshift = (length/2)*cos(atan(slope));
    
    Line *perpendicularline = newPerpendicularLine(inputLine);
	
    operaPoint leftMidPoint(midPoint.getXcoord() - midPointXshift,getYcoord(midPoint.getXcoord() - midPointXshift));
    
    Line *leftline = new Line(perpendicularline->getSlope(),length,width,leftMidPoint);
	
    delete perpendicularline;    
	
    return leftline;
}

Line *Line::newRightLine(const Line &inputLine) const {
    float midPointXshift = (length/2)*cos(atan(slope));
	
    Line *perpendicularline = newPerpendicularLine(inputLine);
    
    operaPoint rightMidPoint(midPoint.getXcoord() + midPointXshift,getYcoord(midPoint.getXcoord() + midPointXshift));
	
    Line *rightline = new Line(perpendicularline->getSlope(),length,width,rightMidPoint);
	
    delete perpendicularline;    
	
    return rightline;
}

operaPoint* Line::newIntersectionPoint(const Line &inputLine) const {
    if(inputLine.getSlope() == getSlope()) {
        throw operaException("Line: error: no intersection point; input line is parallel. ", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);                
    }
    
    float xintersec = (inputLine.getIntercept() - getIntercept())/(getSlope()-inputLine.getSlope());
    float yintersec = inputLine.getYcoord(xintersec);
    
    operaPoint *intersectionPoint = new operaPoint(xintersec,yintersec);
    
    return intersectionPoint;
}


float Line::getIntercept(void) const {
    return intercept;
}

bool Line::pointOnLine(const operaPoint &TestPoint) const {
    float ShiftInXDirection = (length/2)*cos(atan(slope));
	
    float xmin = midPoint.getXcoord() - ShiftInXDirection;
    float xmax = midPoint.getXcoord() + ShiftInXDirection;
	
    if(TestPoint.getYcoord() == getYcoord(TestPoint.getXcoord()) && 
       TestPoint.getXcoord() >= xmin &&
       TestPoint.getXcoord() <= xmax ) {
        return true;
    }
    return false;
}

void Line::setMidPoint(const operaPoint &MidPoint) {
    midPoint.setPoint(MidPoint.getXcoord(), MidPoint.getYcoord());
}

void Line::setMidPoint(const operaPoint *MidPoint) {
    midPoint.setPoint(MidPoint->getXcoord(), MidPoint->getYcoord());
}
bool Line::pointInLineWidth(const operaPoint &TestPoint, const Line &inputLine) const {
    bool PointIsIN = FALSE;
	
    Line *topline = newTopLine(inputLine);
    Line *bottomline = newBottomLine(inputLine);
    Line *leftline = newLeftLine(inputLine);
    Line *rightline = newRightLine(inputLine);
	
    if(TestPoint.getYcoord() < topline->getYcoord(TestPoint.getXcoord()) &&
       TestPoint.getYcoord() >= bottomline->getYcoord(TestPoint.getXcoord()) &&
       TestPoint.getXcoord() >= leftline->getXcoord(TestPoint.getYcoord()) &&
       TestPoint.getXcoord() < rightline->getXcoord(TestPoint.getYcoord()) ) {
        PointIsIN = TRUE;
    }        
	
    delete topline;
    delete bottomline;
    delete leftline;
    delete rightline;
    
    return PointIsIN;
}

float Line::getLineYWidth(void) const {
    float ShiftInYDirection = width/(2*cos(atan(slope)));
    
    return 2*ShiftInYDirection;
}
