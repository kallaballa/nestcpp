#define _USE_MATH_DEFINES
#include "geometry_util.hpp"


namespace GeometryUtil {
bool _almostEqual(dim_t a, dim_t b, dim_t tolerance) {
	return fabs(a - b) < tolerance;
}

bool _withinDistance(const Point& p1, const Point& p2, dim_t distance){
	dim_t dx = p1.x-p2.x;
	dim_t dy = p1.y-p2.y;
	return ((dx*dx + dy*dy) < distance*distance);
}

dim_t _degreesToRadians(dim_t angle){
	return angle*(M_PI/180);
}

dim_t _radiansToDegrees(dim_t angle){
	return angle*(180/M_PI);
}

Point _normalizeVector(const Point& v){
	if(_almostEqual(v.x*v.x + v.y*v.y, 1)){
		return v; // given vector was already a unit vector
	}
	dim_t len = sqrt(v.x*v.x + v.y*v.y);
	dim_t inverse = 1/len;

	return { v.x*inverse,v.y*inverse };
}


bool _onSegment(const Point& A,const Point& B, const Point& p){

	// vertical line
	if(_almostEqual(A.x, B.x) && _almostEqual(p.x, A.x)){
		if(!_almostEqual(p.y, B.y) && !_almostEqual(p.y, A.y) && p.y < std::max(B.y, A.y) && p.y > std::min(B.y, A.y)){
			return true;
		}
		else{
			return false;
		}
	}

	// horizontal line
	if(_almostEqual(A.y, B.y) && _almostEqual(p.y, A.y)){
		if(!_almostEqual(p.x, B.x) && !_almostEqual(p.x, A.x) && p.x < std::max(B.x, A.x) && p.x > std::min(B.x, A.x)){
			return true;
		}
		else{
			return false;
		}
	}

	//range check
	if((p.x < A.x && p.x < B.x) || (p.x > A.x && p.x > B.x) || (p.y < A.y && p.y < B.y) || (p.y > A.y && p.y > B.y)){
		return false;
	}


	// exclude end points
	if((_almostEqual(p.x, A.x) && _almostEqual(p.y, A.y)) || (_almostEqual(p.x, B.x) && _almostEqual(p.y, B.y))){
		return false;
	}

	dim_t cross = (p.y - A.y) * (B.x - A.x) - (p.x - A.x) * (B.y - A.y);

	if(fabs(cross) > FLOAT_TOL){
		return false;
	}

	dim_t dot = (p.x - A.x) * (B.x - A.x) + (p.y - A.y)*(B.y - A.y);



	if(dot < 0 || _almostEqual(dot, 0)){
		return false;
	}

	dim_t len2 = (B.x - A.x)*(B.x - A.x) + (B.y - A.y)*(B.y - A.y);

	if(dot > len2 || _almostEqual(dot, len2)){
		return false;
	}

	return true;
}

Point _lineIntersect(const Point& A,const Point& B,const Point& E,const Point& F, bool infinite){
	dim_t a1, a2, b1, b2, c1, c2, x, y;

	a1= B.y-A.y;
	b1= A.x-B.x;
	c1= B.x*A.y - A.x*B.y;
	a2= F.y-E.y;
	b2= E.x-F.x;
	c2= F.x*E.y - E.x*F.y;

	dim_t denom=a1*b2 - a2*b1;

	x = (b1*c2 - b2*c1)/denom;
	y = (a2*c1 - a1*c2)/denom;

	if(!std::isfinite(x) || !std::isfinite(y)){
		return INVALID_POINT;
	}

	// lines are colinear
	/*var crossABE = (E.y - A.y) * (B.x - A.x) - (E.x - A.x) * (B.y - A.y);
	var crossABF = (F.y - A.y) * (B.x - A.x) - (F.x - A.x) * (B.y - A.y);
	if(_almostEqual(crossABE,0) && _almostEqual(crossABF,0)){
		return null;
	}*/

	if(!infinite){
		// coincident points do not count as intersecting
		if (fabs(A.x-B.x) > FLOAT_TOL && (( A.x < B.x ) ? x < A.x || x > B.x : x > A.x || x < B.x )) return INVALID_POINT;
		if (fabs(A.y-B.y) > FLOAT_TOL && (( A.y < B.y ) ? y < A.y || y > B.y : y > A.y || y < B.y )) return INVALID_POINT;

		if (fabs(E.x-F.x) > FLOAT_TOL && (( E.x < F.x ) ? x < E.x || x > F.x : x > E.x || x < F.x )) return INVALID_POINT;
		if (fabs(E.y-F.y) > FLOAT_TOL && (( E.y < F.y ) ? y < E.y || y > F.y : y > E.y || y < F.y )) return INVALID_POINT;
	}

	return {x, y};
}

namespace QuadradicBezier {

bool isFlat(const Point& p1, const Point& p2, const Point& c1, dim_t tol) {
	tol = 4 * tol * tol;

	dim_t ux = 2 * c1.x - p1.x - p2.x;
	ux *= ux;

	dim_t uy = 2 * c1.y - p1.y - p2.y;
	uy *= uy;

	return (ux + uy <= tol);
}

std::pair<Segment, Segment> subdivide(const Point& p1, const Point& p2,
		const Point& c1, dim_t t) {
	Point mid1 = { p1.x + (c1.x - p1.x) * t, p1.y + (c1.y - p1.y) * t };
	Point mid2 = { c1.x + (p2.x - c1.x) * t, c1.y + (p2.y - c1.y) * t };
	Point mid3 = { mid1.x + (mid2.x - mid1.x) * t, mid1.y + (mid2.y - mid1.y) * t };

	Segment seg1 = { p1, mid3, mid1 };
	Segment seg2 = { mid3, p2, mid2 };

	return {seg1, seg2};
}

list<Point> linearize(const Point& p1, const Point& p2, const Point& c1, dim_t tol) {
	list<Point> finished = { p1 }; // list of points to return
	list<Segment> todo = { { p1, p2, c1 } }; // list of Beziers to divide

	// recursion could stack overflow, loop instead
	while (todo.size() > 0) {
		Segment& segment = todo.front();

		if (isFlat(segment.p1, segment.p2, segment.c1, tol)) { // reached subdivision limit
			finished.push_back( { segment.p2.x, segment.p2.y });
			todo.erase(todo.begin());
		} else {
			auto divided = subdivide(segment.p1, segment.p2, segment.c1, 0.5);
			//AMIR, TESTME: does this really splice correctly?
			todo.erase(todo.begin());
			todo.push_front(divided.second);
			todo.push_front(divided.first);
		}
	}
	return finished;
}
} //QuadraticBezier

namespace CubicBezier {

bool isFlat(const Point& p1, const Point& p2, const Point& c1, const Point& c2,	dim_t tol) {
	tol = 16 * tol * tol;

	dim_t ux = 3 * c1.x - 2 * p1.x - p2.x;
	ux *= ux;

	dim_t uy = 3 * c1.y - 2 * p1.y - p2.y;
	uy *= uy;

	dim_t vx = 3 * c2.x - 2 * p2.x - p1.x;
	vx *= vx;

	dim_t vy = 3 * c2.y - 2 * p2.y - p1.y;
	vy *= vy;

	if (ux < vx) {
		ux = vx;
	}
	if (uy < vy) {
		uy = vy;
	}

	return (ux + uy <= tol);
}

std::pair<Segment,Segment> subdivide(const Point& p1, const Point& p2, const Point& c1, const Point& c2, const dim_t& t){
	Point mid1 = {
		p1.x+(c1.x-p1.x)*t,
		p1.y+(c1.y-p1.y)*t
	};

	Point mid2 = {
		c2.x+(p2.x-c2.x)*t,
		c2.y+(p2.y-c2.y)*t
	};

	Point mid3 = {
		c1.x+(c2.x-c1.x)*t,
		c1.y+(c2.y-c1.y)*t
	};

	Point mida = {
		mid1.x+(mid3.x-mid1.x)*t,
		mid1.y+(mid3.y-mid1.y)*t
	};

	Point midb = {
		mid3.x+(mid2.x-mid3.x)*t,
		mid3.y+(mid2.y-mid3.y)*t
	};

	Point midx = {
		mida.x+(midb.x-mida.x)*t,
		mida.y+(midb.y-mida.y)*t
	};

	Segment seg1 = {p1, midx, mid1, mida};
	Segment seg2 = {midx, p2, midb, mid2};

	return {seg1, seg2};
}

list<Point> linearize(const Point& p1, const Point& p2, const Point& c1, const Point& c2, dim_t tol){
	list<Point> finished = {p1}; // list of points to return
	list<Segment> todo = {{p1, p2, c1, c2}}; // list of Beziers to divide

	// recursion could stack overflow, loop instead

	while(todo.size() > 0){
		Segment segment = todo.front();

		if(isFlat(segment.p1, segment.p2, segment.c1, segment.c2, tol)){ // reached subdivision limit
			finished.push_back({segment.p2.x, segment.p2.y});
			todo.erase(todo.begin());
		}
		else{
			auto divided = subdivide(segment.p1, segment.p2, segment.c1, segment.c2, 0.5);
			//AMIR, TESTME: does this really splice correctly?
			todo.erase(todo.begin());
			todo.push_front(divided.second);
			todo.push_front(divided.first);
		}
	}
	return finished;
}
} //CubeBezier
namespace Arc {

SvgArc centerToSvg(const Point& center, const dim_t& rx, const dim_t& ry, dim_t theta1, const dim_t& extent, const dim_t& angleDegrees){
	dim_t theta2 = theta1 + extent;

	theta1 = _degreesToRadians(theta1);
	theta2 = _degreesToRadians(theta2);
	dim_t angle = _degreesToRadians(angleDegrees);

	dim_t dcos = cos(angle);
	dim_t dsin = sin(angle);

	dim_t t1cos = cos(theta1);
	dim_t t1sin = sin(theta1);

	dim_t t2cos = cos(theta2);
	dim_t t2sin = sin(theta2);

	dim_t x0 = center.x + dcos * rx * t1cos +	(-dsin) * ry * t1sin;
	dim_t y0 = center.y + dsin * rx * t1cos +	dcos * ry * t1sin;

	dim_t x1 = center.x + dcos * rx * t2cos +	(-dsin) * ry * t2sin;
	dim_t y1 = center.y + dsin * rx * t2cos +	dcos * ry * t2sin;

	int8_t largearc = (extent > 180) ? 1 : 0;
	int8_t sweep = (extent > 0) ? 1 : 0;

	return {
		{x0, y0},
		{x1, y1},
		rx,
		ry,
		angle,
		largearc,
		sweep
	};
}

CenterArc svgToCenter(const Point& p1, const Point& p2, dim_t rx,
		dim_t ry, const dim_t& angleDegrees, const int8_t& largearc,
		const int8_t& sweep) {
	Point mid = { 0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y) };
	Point diff = { 0.5 * (p2.x - p1.x), 0.5 * (p2.y - p1.y) };

	dim_t angle = _degreesToRadians(((long)angleDegrees) % 360);

	dim_t c = cos(angle);
	dim_t s = sin(angle);

	dim_t x1 = c * diff.x + s * diff.y;
	dim_t y1 = -s * diff.x + c * diff.y;

	rx = fabs(rx);
	ry = fabs(ry);
	dim_t Prx = rx * rx;
	dim_t Pry = ry * ry;
	dim_t Px1 = x1 * x1;
	dim_t Py1 = y1 * y1;

	dim_t radiiCheck = Px1 / Prx + Py1 / Pry;
	dim_t radiiSqrt = sqrt(radiiCheck);
	if (radiiCheck > 1) {
		rx = radiiSqrt * rx;
		ry = radiiSqrt * ry;
		Prx = rx * rx;
		Pry = ry * ry;
	}

	int8_t sign = (largearc != sweep) ? -1 : 1;
	dim_t sq = ((Prx * Pry) - (Prx * Py1) - (Pry * Px1))
			/ ((Prx * Py1) + (Pry * Px1));

	sq = (sq < 0) ? 0 : sq;

	dim_t coef = sign * sqrt(sq);
	dim_t cx1 = coef * ((rx * y1) / ry);
	dim_t cy1 = coef * -((ry * x1) / rx);

	dim_t cx = mid.x + (c * cx1 - s * cy1);
	dim_t cy = mid.y + (s * cx1 + c * cy1);

	dim_t ux = (x1 - cx1) / rx;
	dim_t uy = (y1 - cy1) / ry;
	dim_t vx = (-x1 - cx1) / rx;
	dim_t vy = (-y1 - cy1) / ry;
	dim_t n = sqrt((ux * ux) + (uy * uy));
	dim_t p = ux;
	sign = (uy < 0) ? -1 : 1;

	dim_t theta = sign * acos(p / n);
	theta = _radiansToDegrees(theta);

	n = sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy));
	p = ux * vx + uy * vy;
	sign = ((ux * vy - uy * vx) < 0) ? -1 : 1;
	dim_t delta = sign * acos(p / n);
	delta = _radiansToDegrees(delta);

	if (sweep == 1 && delta > 0) {
		delta -= 360;
	} else if (sweep == 0 && delta < 0) {
		delta += 360;
	}

	delta = ((long)delta) % 360;
	theta = ((long)theta) % 360;

	return {
		{	cx, cy},
		rx,
		ry,
		theta,
		delta,
		angleDegrees
	};
}

list<Point> linearize(const Point& p1, const Point& p2, dim_t rx, dim_t ry,
		dim_t angle, int8_t largearc, int8_t sweep, dim_t tol) {
	list<Point> finished = { p2 }; // list of points to return
	CenterArc arc = svgToCenter(p1, p2, rx, ry, angle, largearc, sweep);
	list<CenterArc> todo = { arc }; // list of arcs to divide

	// recursion could stack overflow, loop instead
	while (todo.size() > 0) {
		arc = todo.front();

		SvgArc fullarc = centerToSvg(arc.center, arc.rx, arc.ry, arc.theta,
				arc.extent, arc.angle);
		SvgArc subarc = centerToSvg(arc.center, arc.rx, arc.ry, arc.theta,
				0.5 * arc.extent, arc.angle);
		Point& arcmid = subarc.p2;

		Point mid = { 0.5 * (fullarc.p1.x + fullarc.p2.x), 0.5
				* (fullarc.p1.y + fullarc.p2.y) };

		// compare midpoint of line with midpoint of arc
		// this is not 100% accurate, but should be a good heuristic for flatness in most cases
		if (_withinDistance(mid, arcmid, tol)) {
			finished.push_front(fullarc.p2);
			todo.erase(todo.begin());
		} else {
			CenterArc arc1 = { arc.center, arc.rx, arc.ry, arc.theta, 0.5 * arc.extent, arc.angle };
			CenterArc arc2 = { arc.center, arc.rx, arc.ry, arc.theta + 0.5 * arc.extent, 0.5 * arc.extent, arc.angle };
			todo.erase(todo.begin());
			todo.push_front(arc2);
			todo.push_front(arc1);
		}
	}
	return finished;
}
} //Arc

Rect getPolygonBounds(const Polygon& polygon){
	if(polygon.size() < 3){
		return INVALID_RECT;
	}

	dim_t xmin = polygon[0].x;
	dim_t xmax = polygon[0].x;
	dim_t ymin = polygon[0].y;
	dim_t ymax = polygon[0].y;

	for(size_t i=1; i<polygon.size(); i++){
		if(polygon[i].x > xmax){
			xmax = polygon[i].x;
		}
		else if(polygon[i].x < xmin){
			xmin = polygon[i].x;
		}

		if(polygon[i].y > ymax){
			ymax = polygon[i].y;
		}
		else if(polygon[i].y < ymin){
			ymin = polygon[i].y;
		}
	}

	return {
		xmin,
		ymin,
		xmax-xmin,
		ymax-ymin
	};
}

PointInPolygonResult pointInPolygon(const Point& point, const Polygon& polygon){
	if(polygon.size() < 3){
		return INVALID;
	}

	bool inside = false;
	dim_t offsetx = polygon.offsetx;
	dim_t offsety = polygon.offsety;

	for (size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j=i++) {
		dim_t xi = polygon[i].x + offsetx;
		dim_t yi = polygon[i].y + offsety;
		dim_t xj = polygon[j].x + offsetx;
		dim_t yj = polygon[j].y + offsety;

		if(_almostEqual(xi, point.x) && _almostEqual(yi, point.y)){
			return INVALID; // no result
		}

		if(_onSegment({xi, yi}, {xj, yj}, point)){
			return INVALID; // exactly on the segment
		}

		if(_almostEqual(xi, xj) && _almostEqual(yi, yj)){ // ignore very small lines
			continue;
		}

		bool intersect = ((yi > point.y) != (yj > point.y)) && (point.x < (xj - xi) * (point.y - yi) / (yj - yi) + xi);
		if (intersect) inside = !inside;
	}

	return inside ? INSIDE : OUTSIDE;
}

dim_t polygonArea(const Polygon& polygon){
	dim_t area = 0;
	size_t i, j;
	for (i=0, j=polygon.size()-1; i<polygon.size(); j=i++){
		area += (polygon[j].x+polygon[i].x) * (polygon[j].y-polygon[i].y);
	}
	return 0.5*area;
}

bool intersect(const Polygon& A, const Polygon& B){
	dim_t Aoffsetx = A.offsetx;
	dim_t Aoffsety = A.offsety;

	dim_t Boffsetx = B.offsetx;
	dim_t Boffsety = B.offsety;


	for(size_t i=0; i<A.size()-1; i++){
		for(size_t j=0; j<B.size()-1; j++){
			Point a1 = {A[i].x+Aoffsetx, A[i].y+Aoffsety};
			Point a2 = {A[i+1].x+Aoffsetx, A[i+1].y+Aoffsety};
			Point b1 = {B[j].x+Boffsetx, B[j].y+Boffsety};
			Point b2 = {B[j+1].x+Boffsetx, B[j+1].y+Boffsety};

			size_t prevbindex = (j == 0) ? B.size()-1 : j-1;
			size_t prevaindex = (i == 0) ? A.size()-1 : i-1;
			size_t nextbindex = (j+1 == B.size()-1) ? 0 : j+2;
			size_t nextaindex = (i+1 == A.size()-1) ? 0 : i+2;

			// go even further back if we happen to hit on a loop end point
			if(B[prevbindex] == B[j] || (_almostEqual(B[prevbindex].x, B[j].x) && _almostEqual(B[prevbindex].y, B[j].y))){
				prevbindex = (prevbindex == 0) ? B.size()-1 : prevbindex-1;
			}

			if(A[prevaindex] == A[i] || (_almostEqual(A[prevaindex].x, A[i].x) && _almostEqual(A[prevaindex].y, A[i].y))){
				prevaindex = (prevaindex == 0) ? A.size()-1 : prevaindex-1;
			}

			// go even further forward if we happen to hit on a loop end point
			if(B[nextbindex] == B[j+1] || (_almostEqual(B[nextbindex].x, B[j+1].x) && _almostEqual(B[nextbindex].y, B[j+1].y))){
				nextbindex = (nextbindex == B.size()-1) ? 0 : nextbindex+1;
			}

			if(A[nextaindex] == A[i+1] || (_almostEqual(A[nextaindex].x, A[i+1].x) && _almostEqual(A[nextaindex].y, A[i+1].y))){
				nextaindex = (nextaindex == A.size()-1) ? 0 : nextaindex+1;
			}

			Point a0 = {A[prevaindex].x + Aoffsetx, A[prevaindex].y + Aoffsety};
			Point b0 = {B[prevbindex].x + Boffsetx, B[prevbindex].y + Boffsety};

			Point a3 = {A[nextaindex].x + Aoffsetx, A[nextaindex].y + Aoffsety};
			Point b3 = {B[nextbindex].x + Boffsetx, B[nextbindex].y + Boffsety};

			if(_onSegment(a1,a2,b1) || (_almostEqual(a1.x, b1.x) && _almostEqual(a1.y, b1.y))){
				// if a point is on a segment, it could intersect or it could not. Check via the neighboring points
				PointInPolygonResult b0in = pointInPolygon(b0, A);
				PointInPolygonResult b2in = pointInPolygon(b2, A);
				//AMIR: TESTME is this comparison correct in terms of handling INVALID?
				if((b0in == INSIDE && b2in == OUTSIDE) || (b0in == OUTSIDE && b2in == INSIDE)){
					return true;
				}
				else{
					continue;
				}
			}

			if(_onSegment(a1,a2,b2) || (_almostEqual(a2.x, b2.x) && _almostEqual(a2.y, b2.y))){
				// if a point is on a segment, it could intersect or it could not. Check via the neighboring points
				PointInPolygonResult b1in = pointInPolygon(b1, A);
				PointInPolygonResult b3in = pointInPolygon(b3, A);
				//AMIR: TESTME is this comparison correct in terms of handling INVALID?
				if((b1in == INSIDE && b3in == OUTSIDE) || (b1in == OUTSIDE && b3in == INSIDE)){
					return true;
				}
				else{
					continue;
				}
			}

			if(_onSegment(b1,b2,a1) || (_almostEqual(a1.x, b2.x) && _almostEqual(a1.y, b2.y))){
				// if a point is on a segment, it could intersect or it could not. Check via the neighboring points
				PointInPolygonResult a0in = pointInPolygon(a0, B);
				PointInPolygonResult a2in = pointInPolygon(a2, B);
				//AMIR: TESTME is this comparison correct in terms of handling INVALID?
				if((a0in == INSIDE && a2in == OUTSIDE) || (a0in == OUTSIDE && a2in == INSIDE)){
					return true;
				}
				else{
					continue;
				}
			}

			if(_onSegment(b1,b2,a2) || (_almostEqual(a2.x, b1.x) && _almostEqual(a2.y, b1.y))){
				// if a point is on a segment, it could intersect or it could not. Check via the neighboring points
				PointInPolygonResult a1in = pointInPolygon(a1, B);
				PointInPolygonResult a3in = pointInPolygon(a3, B);
				//AMIR: TESTME is this comparison correct in terms of handling INVALID?
				if((a1in == INSIDE && a3in == OUTSIDE) || (a1in == OUTSIDE && a3in == INSIDE)){
					return true;
				}
				else{
					continue;
				}
			}

			Point p = _lineIntersect(b1, b2, a1, a2);

			if(p != INVALID_POINT){
				return true;
			}
		}
	}

	return false;
}


std::vector<Point> polygonEdge(const Polygon& polygon, Point normal){
	if(polygon.size() < 3){
		return {};
	}

	normal = _normalizeVector(normal);

	Point direction = {
		-normal.y,
		normal.x
	};

	// find the max and min points, they will be the endpoints of our edge
	dim_t min = DIM_MAX;
	dim_t max = DIM_MIN;

	std::vector<dim_t> dotproduct;

	for(size_t i=0; i<polygon.size(); i++){
		dim_t dot = polygon[i].x*direction.x + polygon[i].y*direction.y;
		dotproduct.push_back(dot);
		if(min == DIM_MAX || dot < min){
			min = dot;
		}
		if(max == DIM_MIN || dot > max){
			max = dot;
		}
	}

	// there may be multiple vertices with min/max values. In which case we choose the one that is normal-most (eg. left most)
	size_t indexmin = 0;
	size_t indexmax = 0;

	dim_t normalmin = DIM_MAX;
	dim_t normalmax = DIM_MIN;

	for(size_t i=0; i<polygon.size(); i++){
		if(_almostEqual(dotproduct[i] , min)){
			dim_t dot = polygon[i].x*normal.x + polygon[i].y*normal.y;
			if(normalmin == DIM_MAX || dot > normalmin){
				normalmin = dot;
				indexmin = i;
			}
		}
		else if(_almostEqual(dotproduct[i] , max)){
			dim_t dot = polygon[i].x*normal.x + polygon[i].y*normal.y;
			if(normalmax == DIM_MIN || dot > normalmax){
				normalmax = dot;
				indexmax = i;
			}
		}
	}

	// now we have two edges bound by min and max points, figure out which edge faces our direction vector

	size_t indexleft = indexmin-1;
	size_t indexright = indexmin+1;

	if(indexleft < 0){
		indexleft = polygon.size()-1;
	}
	if(indexright >= polygon.size()){
		indexright = 0;
	}

	const Point& minvertex = polygon[indexmin];
	const Point& left = polygon[indexleft];
	const Point& right = polygon[indexright];

	Point leftvector = {
		left.x - minvertex.x,
		left.y - minvertex.y
	};

	Point rightvector = {
		right.x - minvertex.x,
		right.y - minvertex.y
	};

	dim_t dotleft = leftvector.x*direction.x + leftvector.y*direction.y;
	dim_t dotright = rightvector.x*direction.x + rightvector.y*direction.y;

	// -1 = left, 1 = right
	int8_t scandirection = -1;

	if(_almostEqual(dotleft, 0)){
		scandirection = 1;
	}
	else if(_almostEqual(dotright, 0)){
		scandirection = -1;
	}
	else{
		dim_t normaldotleft;
		dim_t normaldotright;

		if(_almostEqual(dotleft, dotright)){
			// the points line up exactly along the normal vector
			normaldotleft = leftvector.x*normal.x + leftvector.y*normal.y;
			normaldotright = rightvector.x*normal.x + rightvector.y*normal.y;
		}
		else if(dotleft < dotright){
			// normalize right vertex so normal projection can be directly compared
			normaldotleft = leftvector.x*normal.x + leftvector.y*normal.y;
			normaldotright = (rightvector.x*normal.x + rightvector.y*normal.y)*(dotleft/dotright);
		}
		else{
			// normalize left vertex so normal projection can be directly compared
			normaldotleft = leftvector.x*normal.x + leftvector.y*normal.y * (dotright/dotleft);
			normaldotright = rightvector.x*normal.x + rightvector.y*normal.y;
		}

		if(normaldotleft > normaldotright){
			scandirection = -1;
		}
		else{
			// technically they could be equal, (ie. the segments bound by left and right points are incident)
			// in which case we'll have to climb up the chain until lines are no longer incident
			// for now we'll just not handle it and assume people aren't giving us garbage input..
			scandirection = 1;
		}
	}

	// connect all points between indexmin and indexmax along the scan direction
	std::vector<Point> edge;
	size_t count = 0;
	size_t i=indexmin;
	while(count < polygon.size()){
		if(i >= polygon.size()){
			i=0;
		}
		else if(i < 0){
			i=polygon.size()-1;
		}

		edge.push_back(polygon[i]);

		if(i == indexmax){
			break;
		}
		i += scandirection;
		count++;
	}

	return edge;
}

dim_t pointLineDistance(const Point& p, const Point& s1, const Point& s2, Point normal, bool s1inclusive, bool s2inclusive){

	normal = _normalizeVector(normal);

	Point dir = {
		normal.y,
		-normal.x
	};

	dim_t pdot = p.x*dir.x + p.y*dir.y;
	dim_t s1dot = s1.x*dir.x + s1.y*dir.y;
	dim_t s2dot = s2.x*dir.x + s2.y*dir.y;

	dim_t pdotnorm = p.x*normal.x + p.y*normal.y;
	dim_t s1dotnorm = s1.x*normal.x + s1.y*normal.y;
	dim_t s2dotnorm = s2.x*normal.x + s2.y*normal.y;


	// point is exactly along the edge in the normal direction
	if(_almostEqual(pdot, s1dot) && _almostEqual(pdot, s2dot)){
		// point lies on an endpoint
		if(_almostEqual(pdotnorm, s1dotnorm) ){
			return DIM_MAX;
		}

		if(_almostEqual(pdotnorm, s2dotnorm)){
			return DIM_MAX;
		}

		// point is outside both endpoints
		if (pdotnorm>s1dotnorm && pdotnorm>s2dotnorm){
			return std::min(pdotnorm - s1dotnorm, pdotnorm - s2dotnorm);
		}
		if (pdotnorm<s1dotnorm && pdotnorm<s2dotnorm){
			return -std::min(s1dotnorm-pdotnorm, s2dotnorm-pdotnorm);
		}

		// point lies between endpoints
		dim_t diff1 = pdotnorm - s1dotnorm;
		dim_t diff2 = pdotnorm - s2dotnorm;
		if(diff1 > 0){
			return diff1;
		}
		else{
			return diff2;
		}
	}
	// point
	else if(_almostEqual(pdot, s1dot)){
		if(s1inclusive){
			return pdotnorm-s1dotnorm;
		}
		else{
			return DIM_MAX;
		}
	}
	else if(_almostEqual(pdot, s2dot)){
		if(s2inclusive){
			return pdotnorm-s2dotnorm;
		}
		else{
			return DIM_MAX;
		}
	}
	else if ((pdot<s1dot && pdot<s2dot) || (pdot>s1dot && pdot>s2dot)){
		return DIM_MAX; // point doesn't collide with segment
	}

	return (pdotnorm - s1dotnorm + (s1dotnorm - s2dotnorm)*(s1dot - pdot)/(s1dot - s2dot));
}

dim_t pointDistance(const Point& p, const Point& s1, const Point& s2,
		Point normal, bool infinite) {
	normal = _normalizeVector(normal);

	Point dir = { normal.y, -normal.x };

	dim_t pdot = p.x * dir.x + p.y * dir.y;
	dim_t s1dot = s1.x * dir.x + s1.y * dir.y;
	dim_t s2dot = s2.x * dir.x + s2.y * dir.y;

	dim_t pdotnorm = p.x * normal.x + p.y * normal.y;
	dim_t s1dotnorm = s1.x * normal.x + s1.y * normal.y;
	dim_t s2dotnorm = s2.x * normal.x + s2.y * normal.y;

	if (!infinite) {
		if (((pdot < s1dot || _almostEqual(pdot, s1dot))
				&& (pdot < s2dot || _almostEqual(pdot, s2dot)))
				|| ((pdot > s1dot || _almostEqual(pdot, s1dot))
						&& (pdot > s2dot || _almostEqual(pdot, s2dot)))) {
			return DIM_MAX; // dot doesn't collide with segment, or lies directly on the vertex
		}
		if ((_almostEqual(pdot, s1dot) && _almostEqual(pdot, s2dot))
				&& (pdotnorm > s1dotnorm && pdotnorm > s2dotnorm)) {
			return std::min(pdotnorm - s1dotnorm, pdotnorm - s2dotnorm);
		}
		if ((_almostEqual(pdot, s1dot) && _almostEqual(pdot, s2dot))
				&& (pdotnorm < s1dotnorm && pdotnorm < s2dotnorm)) {
			return -std::min(s1dotnorm - pdotnorm, s2dotnorm - pdotnorm);
		}
	}

	return -(pdotnorm - s1dotnorm
			+ (s1dotnorm - s2dotnorm) * (s1dot - pdot) / (s1dot - s2dot));
}

dim_t segmentDistance(const Point& A, const Point& B, const Point& E,
		const Point& F, const Point& direction) {
	Point normal = { direction.y, -direction.x };

	Point reverse = { -direction.x, -direction.y };

	dim_t dotA = A.x * normal.x + A.y * normal.y;
	dim_t dotB = B.x * normal.x + B.y * normal.y;
	dim_t dotE = E.x * normal.x + E.y * normal.y;
	dim_t dotF = F.x * normal.x + F.y * normal.y;

	dim_t crossA = A.x * direction.x + A.y * direction.y;
	dim_t crossB = B.x * direction.x + B.y * direction.y;
	dim_t crossE = E.x * direction.x + E.y * direction.y;
	dim_t crossF = F.x * direction.x + F.y * direction.y;

	dim_t ABmin = std::min(dotA, dotB);
	dim_t ABmax = std::max(dotA, dotB);

	dim_t EFmax = std::max(dotE, dotF);
	dim_t EFmin = std::min(dotE, dotF);

	// segments that will merely touch at one point
	if (_almostEqual(ABmax, EFmin, FLOAT_TOL)
			|| _almostEqual(ABmin, EFmax, FLOAT_TOL)) {
		return DIM_MAX;
	}
	// segments miss eachother completely
	if (ABmax < EFmin || ABmin > EFmax) {
		return DIM_MAX;
	}

	dim_t overlap;

	if ((ABmax > EFmax && ABmin < EFmin) || (EFmax > ABmax && EFmin < ABmin)) {
		overlap = 1;
	} else {
		dim_t minMax = std::min(ABmax, EFmax);
		dim_t maxMin = std::max(ABmin, EFmin);

		dim_t maxMax = std::max(ABmax, EFmax);
		dim_t minMin = std::min(ABmin, EFmin);

		overlap = (minMax - maxMin) / (maxMax - minMin);
	}

	dim_t crossABE = (E.y - A.y) * (B.x - A.x) - (E.x - A.x) * (B.y - A.y);
	dim_t crossABF = (F.y - A.y) * (B.x - A.x) - (F.x - A.x) * (B.y - A.y);

	// lines are colinear
	if (_almostEqual(crossABE, 0) && _almostEqual(crossABF, 0)) {

		Point ABnorm = { B.y - A.y, A.x - B.x };
		Point EFnorm = { F.y - E.y, E.x - F.x };

		dim_t ABnormlength = sqrt(ABnorm.x * ABnorm.x + ABnorm.y * ABnorm.y);
		ABnorm.x /= ABnormlength;
		ABnorm.y /= ABnormlength;

		dim_t EFnormlength = sqrt(EFnorm.x * EFnorm.x + EFnorm.y * EFnorm.y);
		EFnorm.x /= EFnormlength;
		EFnorm.y /= EFnormlength;

		// segment normals must point in opposite directions
		if (fabs(ABnorm.y * EFnorm.x - ABnorm.x * EFnorm.y) < FLOAT_TOL
				&& ABnorm.y * EFnorm.y + ABnorm.x * EFnorm.x < 0) {
			// normal of AB segment must point in same direction as given direction vector
			dim_t normdot = ABnorm.y * direction.y + ABnorm.x * direction.x;
			// the segments merely slide along eachother
			if (_almostEqual(normdot, 0, FLOAT_TOL)) {
				return DIM_MAX;
			}
			if (normdot < 0) {
				return 0;
			}
		}
		return DIM_MAX;
	}

	std::vector<dim_t> distances;

	// coincident points
	if (_almostEqual(dotA, dotE)) {
		distances.push_back(crossA - crossE);
	} else if (_almostEqual(dotA, dotF)) {
		distances.push_back(crossA - crossF);
	} else if (dotA > EFmin && dotA < EFmax) {
		dim_t d = pointDistance(A, E, F, reverse);
		if (d != DIM_MAX && _almostEqual(d, 0)) { //  A currently touches EF, but AB is moving away from EF
			dim_t dB = pointDistance(B, E, F, reverse, true);
			if (dB < 0 || _almostEqual(dB * overlap, 0)) {
				d = DIM_MAX;
			}
		}
		if (d != DIM_MAX) {
			distances.push_back(d);
		}
	}

	if (_almostEqual(dotB, dotE)) {
		distances.push_back(crossB - crossE);
	} else if (_almostEqual(dotB, dotF)) {
		distances.push_back(crossB - crossF);
	} else if (dotB > EFmin && dotB < EFmax) {
		dim_t d = pointDistance(B, E, F, reverse);

		if (d != DIM_MAX && _almostEqual(d, 0)) { // crossA>crossB A currently touches EF, but AB is moving away from EF
			dim_t dA = pointDistance(A, E, F, reverse, true);
			if (dA < 0 || _almostEqual(dA * overlap, 0)) {
				d = DIM_MAX;
			}
		}
		if (d != DIM_MAX) {
			distances.push_back(d);
		}
	}

	if (dotE > ABmin && dotE < ABmax) {
		dim_t d = pointDistance(E, A, B, direction);
		if (d != DIM_MAX && _almostEqual(d, 0)) { // crossF<crossE A currently touches EF, but AB is moving away from EF
			dim_t dF = pointDistance(F, A, B, direction, true);
			if (dF < 0 || _almostEqual(dF * overlap, 0)) {
				d = DIM_MAX;
			}
		}
		if (d != DIM_MAX) {
			distances.push_back(d);
		}
	}

	if (dotF > ABmin && dotF < ABmax) {
		dim_t d = pointDistance(F, A, B, direction);
		if (d != DIM_MAX && _almostEqual(d, 0)) { // && crossE<crossF A currently touches EF, but AB is moving away from EF
			dim_t dE = pointDistance(E, A, B, direction, true);
			if (dE < 0 || _almostEqual(dE * overlap, 0)) {
				d = DIM_MAX;
			}
		}
		if (d != DIM_MAX) {
			distances.push_back(d);
		}
	}

	if (distances.size() == 0) {
		return DIM_MAX;
	}

	return *std::min_element(distances.begin(), distances.end());
}

dim_t polygonSlideDistance(Polygon A, Polygon B, const Point& direction,
		bool ignoreNegative) {

	Point A1, A2, B1, B2;
	dim_t Aoffsetx, Aoffsety, Boffsetx, Boffsety;

	Aoffsetx = A.offsetx;
	Aoffsety = A.offsety;

	Boffsetx = B.offsetx;
	Boffsety = B.offsety;

	// close the loop for polygons
	if (A[0] != A[A.size() - 1]) {
		A.push_back(A.front());
	}

	if (B[0] != B[B.size() - 1]) {
		B.push_back(B.front());
	}

	Polygon edgeA = A;
	Polygon edgeB = B;

	dim_t distance = DIM_MAX;
	dim_t d;

	Point dir = _normalizeVector(direction);

	for (size_t i = 0; i < edgeB.size() - 1; i++) {
		for (size_t j = 0; j < edgeA.size() - 1; j++) {
			A1 = {edgeA[j].x + Aoffsetx, edgeA[j].y + Aoffsety};
			A2 = {edgeA[j+1].x + Aoffsetx, edgeA[j+1].y + Aoffsety};
			B1 = {edgeB[i].x + Boffsetx, edgeB[i].y + Boffsety};
			B2 = {edgeB[i+1].x + Boffsetx, edgeB[i+1].y + Boffsety};

			if((_almostEqual(A1.x, A2.x) && _almostEqual(A1.y, A2.y)) || (_almostEqual(B1.x, B2.x) && _almostEqual(B1.y, B2.y))) {
				continue; // ignore extremely small lines
			}

			d = segmentDistance(A1, A2, B1, B2, dir);

			if(d != DIM_MAX && (distance == DIM_MAX || d < distance)) {
				if(!ignoreNegative || d > 0 || _almostEqual(d, 0.0)) {
					distance = d;
				}
			}
		}
	}
	return distance;
}

dim_t polygonProjectionDistance(Polygon A, Polygon B, const Point& direction) {
	dim_t Boffsetx = B.offsetx;
	dim_t Boffsety = B.offsety;

	dim_t Aoffsetx = A.offsetx;
	dim_t Aoffsety = A.offsety;

	// close the loop for polygons
	if (A[0] != A[A.size() - 1]) {
		A.push_back(A[0]);
	}

	if (B[0] != B[B.size() - 1]) {
		B.push_back(B[0]);
	}

	Polygon edgeA = A;
	Polygon edgeB = B;

	dim_t distance = DIM_MAX;
	Point p, s1, s2;
	dim_t d;

	for (size_t i = 0; i < edgeB.size(); i++) {
		// the shortest/most negative projection of B onto A
		dim_t minprojection = DIM_MAX;
		Point minp = INVALID_POINT;
		for (size_t j = 0; j < edgeA.size() - 1; j++) {
			p = {edgeB[i].x + Boffsetx, edgeB[i].y + Boffsety};
			s1 = {edgeA[j].x + Aoffsetx, edgeA[j].y + Aoffsety};
			s2 = {edgeA[j+1].x + Aoffsetx, edgeA[j+1].y + Aoffsety};

			if(fabs((s2.y-s1.y) * direction.x - (s2.x-s1.x) * direction.y) < FLOAT_TOL) {
				continue;
			}

			// project point, ignore edge boundaries
			d = pointDistance(p, s1, s2, direction);

			if(d != DIM_MAX && (minprojection == DIM_MAX || d < minprojection)) {
				minprojection = d;
				minp = p;
			}
		}
		if (minprojection != DIM_MAX
				&& (distance == DIM_MAX || minprojection > distance)) {
			distance = minprojection;
		}
	}

	return distance;
}

// returns true if point already exists in the given nfp
bool inNfp(const Point& p, const std::vector<Polygon>& nfp) {
	if (nfp.empty()) {
		return false;
	}

	for (size_t i = 0; i < nfp.size(); i++) {
		for (size_t j = 0; j < nfp[i].size(); j++) {
			if (_almostEqual(p.x, nfp[i][j].x) && _almostEqual(p.y, nfp[i][j].y)) {
				return true;
			}
		}
	}

	return false;
}
Point searchStartPoint(Polygon A, Polygon B, bool inside,
		const vector<Polygon>& NFP) {

	// close the loop for polygons
	if (A[0] != A[A.size() - 1]) {
		A.push_back(A[0]);
	}

	if (B[0] != B[B.size() - 1]) {
		B.push_back(B[0]);
	}

	for (size_t i = 0; i < A.size() - 1; i++) {
		if (!A[i].marked) {
			A[i].marked = true;
			for (size_t j = 0; j < B.size(); j++) {
				B.offsetx = A[i].x - B[j].x;
				B.offsety = A[i].y - B[j].y;

				PointInPolygonResult Binside = INVALID;
				for (size_t k = 0; k < B.size(); k++) {
					PointInPolygonResult inpoly = pointInPolygon( { B[k].x + B.offsetx,
							 B[k].y + B.offsety }, A);
					if (inpoly != INVALID) {
						Binside = inpoly;
						break;
					}
				}

				if (Binside == INVALID) { // A and B are the same
					return INVALID_POINT;
				}

				Point startPoint = { B.offsetx, B.offsety };
				if (((Binside && inside) || (!Binside && !inside)) && !intersect(A, B)
						&& !inNfp(startPoint, NFP)) {
					return startPoint;
				}

				// slide B along vector
				dim_t vx = A[i + 1].x - A[i].x;
				dim_t vy = A[i + 1].y - A[i].y;

				dim_t d1 = polygonProjectionDistance(A, B, {vx, vy });
				dim_t d2 = polygonProjectionDistance(B, A, {-vx, -vy });

				dim_t d = DIM_MAX;

				// todo: clean this up
				if (d1 == DIM_MAX && d2 == DIM_MAX) {
					// nothin
				} else if (d1 == DIM_MAX) {
					d = d2;
				} else if (d2 == DIM_MAX) {
					d = d1;
				} else {
					d = std::min(d1, d2);
				}

				// only slide until no longer negative
				// todo: clean this up
				if (d == DIM_MAX && !_almostEqual(d, 0) && d > 0) {

				} else {
					continue;
				}

				dim_t vd2 = vx * vx + vy * vy;

				if (d * d < vd2 && !_almostEqual(d * d, vd2)) {
					dim_t vd = sqrt(vx * vx + vy * vy);
					vx *= d / vd;
					vy *= d / vd;
				}

				B.offsetx += vx;
				B.offsety += vy;

				for (size_t k = 0; k < B.size(); k++) {
					PointInPolygonResult inpoly = pointInPolygon( { B[k].x + B.offsetx,
							B[k].y + B.offsety }, A);
					if (inpoly == INVALID) {
						Binside = inpoly;
						break;
					}
				}
				startPoint = {B.offsetx, B.offsety};
				if (((Binside && inside) || (!Binside && !inside)) && !intersect(A, B)
						&& !inNfp(startPoint, NFP)) {
					return startPoint;
				}
			}
		}
	}

	return INVALID_POINT;
}

bool isRectangle(const Polygon& poly, dim_t tolerance) {
	Rect bb = getPolygonBounds(poly);

	for (size_t i = 0; i < poly.size(); i++) {
		if (!_almostEqual(poly[i].x, bb.x)
				&& !_almostEqual(poly[i].x, bb.x + bb.width)) {
			return false;
		}
		if (!_almostEqual(poly[i].y, bb.y)
				&& !_almostEqual(poly[i].y, bb.y + bb.height)) {
			return false;
		}
	}

	return true;
}

// returns an interior NFP for the special case where A is a rectangle
vector<Polygon> noFitPolygonRectangle(const Polygon& A, const Polygon& B) {
	dim_t minAx = A[0].x;
	dim_t minAy = A[0].y;
	dim_t maxAx = A[0].x;
	dim_t maxAy = A[0].y;

	for (size_t i = 1; i < A.size(); i++) {
		if (A[i].x < minAx) {
			minAx = A[i].x;
		}
		if (A[i].y < minAy) {
			minAy = A[i].y;
		}
		if (A[i].x > maxAx) {
			maxAx = A[i].x;
		}
		if (A[i].y > maxAy) {
			maxAy = A[i].y;
		}
	}

	dim_t minBx = B[0].x;
	dim_t minBy = B[0].y;
	dim_t maxBx = B[0].x;
	dim_t maxBy = B[0].y;

	for (size_t i = 1; i < B.size(); i++) {
		if (B[i].x < minBx) {
			minBx = B[i].x;
		}
		if (B[i].y < minBy) {
			minBy = B[i].y;
		}
		if (B[i].x > maxBx) {
			maxBx = B[i].x;
		}
		if (B[i].y > maxBy) {
			maxBy = B[i].y;
		}
	}

	if (maxBx - minBx > maxAx - minAx) {
		return {};
	}
	if (maxBy - minBy > maxAy - minAy) {
		return {};
	}

	return { {
			{	minAx-minBx+B[0].x, minAy-minBy+B[0].y},
			{	maxAx-maxBx+B[0].x, minAy-minBy+B[0].y},
			{	maxAx-maxBx+B[0].x, maxAy-maxBy+B[0].y},
			{	minAx-minBx+B[0].x, maxAy-maxBy+B[0].y}
		}};
}

// given a static polygon A and a movable polygon B, compute a no fit polygon by orbiting B about A
	// if the inside flag is set, B is orbited inside of A rather than outside
	// if the searchEdges flag is set, all edges of A are explored for NFPs - multiple
vector<Polygon>	noFitPolygon(Polygon A, Polygon B, bool inside, bool searchEdges){
		if(A.size() < 3 || B.size() < 3){
			return {};
		}

		A.offsetx = 0;
		A.offsety = 0;

		size_t i, j;

		dim_t minA = A[0].y;
		size_t minAindex = 0;

		dim_t maxB = B[0].y;
		size_t maxBindex = 0;

		for(i=1; i<A.size(); i++){
			A[i].marked = false;
			if(A[i].y < minA){
				minA = A[i].y;
				minAindex = i;
			}
		}

		for(i=1; i<B.size(); i++){
			B[i].marked = false;
			if(B[i].y > maxB){
				maxB = B[i].y;
				maxBindex = i;
			}
		}

		Point startpoint;
		if(!inside){
			// shift B such that the bottom-most point of B is at the top-most point of A. This guarantees an initial placement with no intersections
			startpoint = {
				A[minAindex].x-B[maxBindex].x,
				A[minAindex].y-B[maxBindex].y
			};
		}
		else{
			// no reliable heuristic for inside
			startpoint = searchStartPoint(A,B,true);
		}

		vector<Polygon> NFPlist;

		struct EdgeDescriptor {
			int8_t type;
			size_t A;
			size_t B;
		};

		while(startpoint != INVALID_POINT){
			B.offsetx = startpoint.x;
			B.offsety = startpoint.y;


			// maintain a list of touching points/edges
			vector<EdgeDescriptor> touching;

			Point prevvector = INVALID_POINT; // keep track of previous vector
			Polygon NFP;
			NFP.push_back({B[0].x+B.offsetx,B[0].y+B.offsety});

			dim_t referencex = B[0].x+B.offsetx;
			dim_t referencey = B[0].y+B.offsety;
			dim_t startx = referencex;
			dim_t starty = referencey;
			size_t counter = 0;

			while(counter < 10*(A.size() + B.size())){ // sanity check, prevent infinite loop
				touching = {};
				// find touching vertices/edges
				for(i=0; i<A.size(); i++){
					size_t nexti = (i==A.size()-1) ? 0 : i+1;
					for(j=0; j<B.size(); j++){
						size_t nextj = (j==B.size()-1) ? 0 : j+1;
						if(_almostEqual(A[i].x, B[j].x+B.offsetx) && _almostEqual(A[i].y, B[j].y+B.offsety)){
							touching.push_back({ 0, i, j });
						}
						else if(_onSegment(A[i],A[nexti],{B[j].x+B.offsetx, B[j].y + B.offsety})){
							touching.push_back({	1, nexti, j });
						}
						else if(_onSegment({B[j].x+B.offsetx, B[j].y + B.offsety},{B[nextj].x+B.offsetx, B[nextj].y + B.offsety},A[i])){
							touching.push_back({	2, i, nextj });
						}
					}
				}


				struct TransVector {
					dim_t x;
					dim_t y;
					Point start;
					Point end;

					bool operator==(const TransVector& other) const {
						return this->x == other.x && this->y == other.y && this->start == other.start && this->end == other.end;
					}

					bool operator!=(const TransVector& other) const {
						return !this->operator==(other);
					}
				};

				const TransVector INVALID_TRANSVECTOR = {DIM_MAX, DIM_MAX, INVALID_POINT, INVALID_POINT};
				// generate translation vectors from touching vertices/edges
				vector<TransVector> vectors;
				for(i=0; i<touching.size(); i++){
					Point vertexA = A[touching[i].A];
					vertexA.marked = true;

					// adjacent A vertices
					size_t prevAindex = touching[i].A-1;
					size_t nextAindex = touching[i].A+1;

					prevAindex = (prevAindex < 0) ? A.size()-1 : prevAindex; // loop
					nextAindex = (nextAindex >= A.size()) ? 0 : nextAindex; // loop

					Point prevA = A[prevAindex];
					Point nextA = A[nextAindex];

					// adjacent B vertices
					Point vertexB = B[touching[i].B];

					size_t prevBindex = touching[i].B-1;
					size_t nextBindex = touching[i].B+1;

					prevBindex = (prevBindex < 0) ? B.size()-1 : prevBindex; // loop
					nextBindex = (nextBindex >= B.size()) ? 0 : nextBindex; // loop

					Point prevB = B[prevBindex];
					Point nextB = B[nextBindex];

					if(touching[i].type == 0){
						TransVector vA1 = {
							prevA.x-vertexA.x,
							prevA.y-vertexA.y,
							vertexA,
							prevA
						};

						TransVector vA2 = {
							nextA.x-vertexA.x,
							nextA.y-vertexA.y,
							vertexA,
							nextA
						};

						// B vectors need to be inverted
						TransVector vB1 = {
							vertexB.x-prevB.x,
							vertexB.y-prevB.y,
							prevB,
							vertexB
						};

						TransVector vB2 = {
							vertexB.x-nextB.x,
							vertexB.y-nextB.y,
							nextB,
							vertexB
						};

						vectors.push_back(vA1);
						vectors.push_back(vA2);
						vectors.push_back(vB1);
						vectors.push_back(vB2);
					}
					else if(touching[i].type == 1){
						vectors.push_back({
							vertexA.x-(vertexB.x+B.offsetx),
							vertexA.y-(vertexB.y+B.offsety),
							prevA,
							vertexA
						});

						vectors.push_back({
							prevA.x-(vertexB.x+B.offsetx),
							prevA.y-(vertexB.y+B.offsety),
							vertexA,
							prevA
						});
					}
					else if(touching[i].type == 2){
						vectors.push_back({
							vertexA.x-(vertexB.x+B.offsetx),
							vertexA.y-(vertexB.y+B.offsety),
							prevB,
							vertexB
						});

						vectors.push_back({
							vertexA.x-(prevB.x+B.offsetx),
							vertexA.y-(prevB.y+B.offsety),
							vertexB,
							prevB
						});
					}
				}

				// todo: there should be a faster way to reject vectors that will cause immediate intersection. For now just check them all

				TransVector translate = INVALID_TRANSVECTOR;
				dim_t maxd = 0;

				for(i=0; i<vectors.size(); i++){
					if(vectors[i].x == 0 && vectors[i].y == 0){
						continue;
					}

					// if this vector points us back to where we came from, ignore it.
					// ie cross product = 0, dot product < 0
					if(prevvector != INVALID_POINT && vectors[i].y * prevvector.y + vectors[i].x * prevvector.x < 0){
						// compare magnitude with unit vectors
						dim_t vectorlength = sqrt(vectors[i].x*vectors[i].x+vectors[i].y*vectors[i].y);
						Point unitv = {vectors[i].x/vectorlength, vectors[i].y/vectorlength};

						dim_t prevlength = sqrt(prevvector.x*prevvector.x+prevvector.y*prevvector.y);
						Point prevunit = {prevvector.x/prevlength, prevvector.y/prevlength};

						// we need to scale down to unit vectors to normalize vector length. Could also just do a tan here
						if(fabs(unitv.y * prevunit.x - unitv.x * prevunit.y) < 0.0001){
							continue;
						}
					}

					Point pv = {vectors[i].x, vectors[i].y};
					dim_t d = polygonSlideDistance(A, B, pv, true);
					dim_t vecd2 = vectors[i].x*vectors[i].x + vectors[i].y*vectors[i].y;

					if(d == DIM_MAX || d*d > vecd2){
						dim_t vecd = sqrt(vectors[i].x*vectors[i].x + vectors[i].y*vectors[i].y);
						d = vecd;
					}

					if(d != DIM_MAX && d > maxd){
						maxd = d;
						translate = vectors[i];
					}
				}


				if(translate == INVALID_TRANSVECTOR || _almostEqual(maxd, 0)){
					// didn't close the loop, something went wrong here
					NFP = {};
					break;
				}

				translate.start.marked = true;
				translate.end.marked = true;

				prevvector = Point(translate.x, translate.y);

				// trim
				dim_t vlength2 = translate.x*translate.x + translate.y*translate.y;
				if(maxd*maxd < vlength2 && !_almostEqual(maxd*maxd, vlength2)){
					dim_t scale = sqrt((maxd*maxd)/vlength2);
					translate.x *= scale;
					translate.y *= scale;
				}

				referencex += translate.x;
				referencey += translate.y;

				if(_almostEqual(referencex, startx) && _almostEqual(referencey, starty)){
					// we've made a full loop
					break;
				}

				// if A and B start on a touching horizontal line, the end point may not be the start point
				bool looped = false;
				if(!NFP.empty()){
					for(i=0; i<NFP.size()-1; i++){
						if(_almostEqual(referencex, NFP[i].x) && _almostEqual(referencey, NFP[i].y)){
							looped = true;
						}
					}
				}

				if(looped){
					// we've made a full loop
					break;
				}

				NFP.push_back({
					referencex,
					referencey
				});

				B.offsetx += translate.x;
				B.offsety += translate.y;

				counter++;
			}

			if(!NFP.empty()){
				NFPlist.push_back(NFP);
			}

			if(!searchEdges){
				// only get outer NFP or first inner NFP
				break;
			}

			startpoint = searchStartPoint(A,B,inside,NFPlist);
		}

		return NFPlist;
	}

// given two polygons that touch at at least one point, but do not intersect. Return the outer perimeter of both polygons as a single continuous polygon
// A and B must have the same winding direction
Polygon polygonHull(Polygon A, Polygon B) {
	if (A.size() < 3 || B.size() < 3) {
		return {};
	}

	size_t i, j;

	dim_t Aoffsetx = A.offsetx;
	dim_t Aoffsety = A.offsety;
	dim_t Boffsetx = B.offsetx;
	dim_t Boffsety = B.offsety;

	// start at an extreme point that is guaranteed to be on the final polygon
	dim_t miny = A[0].y;
	//FIXME: AMIR use pointers
	Polygon startPolygon = A;
	size_t startIndex = 0;

	for (i = 0; i < A.size(); i++) {
		if (A[i].y + Aoffsety < miny) {
			miny = A[i].y + Aoffsety;
			startPolygon = A;
			startIndex = i;
		}
	}

	for (i = 0; i < B.size(); i++) {
		if (B[i].y + Boffsety < miny) {
			miny = B[i].y + Boffsety;
			startPolygon = B;
			startIndex = i;
		}
	}

	// for simplicity we'll define polygon A as the starting polygon
	if (startPolygon == B) {
		B = A;
		A = startPolygon;
		Aoffsetx = A.offsetx || 0;
		Aoffsety = A.offsety || 0;
		Boffsetx = B.offsetx || 0;
		Boffsety = B.offsety || 0;
	}

	list<Point> C;
	size_t current = startIndex;
	size_t intercept1 = std::numeric_limits<size_t>::max();
	size_t intercept2 = std::numeric_limits<size_t>::max();

	// scan forward from the starting point
	for (i = 0; i < A.size() + 1; i++) {
		current = (current == A.size()) ? 0 : current;
		size_t next = (current == A.size() - 1) ? 0 : current + 1;
		bool touching = false;
		for (j = 0; j < B.size(); j++) {
			size_t nextj = (j == B.size() - 1) ? 0 : j + 1;
			if (_almostEqual(A[current].x + Aoffsetx, B[j].x + Boffsetx)
					&& _almostEqual(A[current].y + Aoffsety, B[j].y + Boffsety)) {
				C.push_back( { A[current].x + Aoffsetx, A[current].y + Aoffsety });
				intercept1 = j;
				touching = true;
				break;
			} else if (_onSegment(
					{ A[current].x + Aoffsetx, A[current].y + Aoffsety },
					{ A[next].x + Aoffsetx, A[next].y + Aoffsety },
					{ B[j].x + Boffsetx, B[j].y + Boffsety })) {
				C.push_back( { A[current].x + Aoffsetx, A[current].y + Aoffsety });
				C.push_back( { B[j].x + Boffsetx, B[j].y + Boffsety });
				intercept1 = j;
				touching = true;
				break;
			} else if (_onSegment( { B[j].x + Boffsetx, B[j].y + Boffsety },
					{B[nextj].x + Boffsetx, B[nextj].y + Boffsety },
					{ A[current].x + Aoffsetx, A[current].y + Aoffsety })) {
				C.push_back( { A[current].x + Aoffsetx, A[current].y + Aoffsety });
				C.push_back( { B[nextj].x + Boffsetx, B[nextj].y + Boffsety });
				intercept1 = nextj;
				touching = true;
				break;
			}
		}

		if (touching) {
			break;
		}

		C.push_back( { A[current].x + Aoffsetx, A[current].y + Aoffsety });

		current++;
	}

	// scan backward from the starting point
	current = startIndex - 1;
	for (i = 0; i < A.size() + 1; i++) {
		current = (current < 0) ? A.size() - 1 : current;
		size_t next = (current == 0) ? A.size() - 1 : current - 1;
		bool touching = false;
		for (j = 0; j < B.size(); j++) {
			size_t nextj = (j == B.size() - 1) ? 0 : j + 1;
			if (_almostEqual(A[current].x + Aoffsetx, B[j].x + Boffsetx)
					&& _almostEqual(A[current].y, B[j].y + Boffsety)) {
				C.push_front( { A[current].x + Aoffsetx, A[current].y + Aoffsety });
				intercept2 = j;
				touching = true;
				break;
			} else if (_onSegment(
					{ A[current].x + Aoffsetx, A[current].y + Aoffsety },
					{ A[next].x + Aoffsetx, A[next].y + Aoffsety },
					{ B[j].x + Boffsetx, B[j].y + Boffsety })) {
				C.push_front( { A[current].x + Aoffsetx, A[current].y + Aoffsety });
				C.push_front( { B[j].x + Boffsetx, B[j].y + Boffsety });
				intercept2 = j;
				touching = true;
				break;
			} else if (_onSegment(
					{ B[j].x + Boffsetx, B[j].y + Boffsety },
					{ B[nextj].x + Boffsetx, B[nextj].y + Boffsety },
					{ A[current].x + Aoffsetx, A[current].y + Aoffsety })) {
				C.push_front( { A[current].x + Aoffsetx, A[current].y + Aoffsety });
				intercept2 = j;
				touching = true;
				break;
			}
		}

		if (touching) {
			break;
		}

		C.push_front( { A[current].x + Aoffsetx, A[current].y + Aoffsety });

		current--;
	}

	if (intercept1 == std::numeric_limits<size_t>::max()
			|| intercept2 == std::numeric_limits<size_t>::max()) {
		// polygons not touching?
		return {};
	}

	// the relevant points on B now lie between intercept1 and intercept2
	current = intercept1 + 1;
	for (i = 0; i < B.size(); i++) {
		current = (current == B.size()) ? 0 : current;
		C.push_back( { B[current].x + Boffsetx, B[current].y + Boffsety });

		if (current == intercept2) {
			break;
		}

		current++;
	}

	// dedupe
	//FIXME: AMIR we need to copy the list because list supports push_front and vector supports access by index
	Polygon PC;
	std::copy(std::begin(C), std::end(C), std::back_inserter(PC));

	for (i = 0; i < PC.size(); i++) {
		size_t next = (i == PC.size() - 1) ? 0 : i + 1;
		if (_almostEqual(PC[i].x, PC[next].x)
				&& _almostEqual(PC[i].y, PC[next].y)) {
			PC.erase(PC.begin() + i);
			i--;
		}
	}

	return PC;
}

Polygon rotatePolygon(const Polygon& polygon, dim_t angle){
	Polygon rotated;;
	angle = angle * M_PI / 180;
	for(size_t i=0; i<polygon.size(); i++){
		dim_t x = polygon[i].x;
		dim_t y = polygon[i].y;
		dim_t x1 = x*cos(angle)-y*sin(angle);
		dim_t y1 = x*sin(angle)+y*cos(angle);

		rotated.push_back({x1, y1});
	}
	// reset bounding box
	Rect bounds = getPolygonBounds(rotated);
	rotated.x = bounds.x;
	rotated.y = bounds.y;
	rotated.width = bounds.width;
	rotated.height = bounds.height;

	return rotated;
}
} //GeometryUtil



int main() {

}
