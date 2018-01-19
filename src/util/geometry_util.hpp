#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>
#include <list>
#include <cstdint>
#include <algorithm>
#include <initializer_list>


namespace GeometryUtil {

using std::list;
using std::vector;

typedef double dim_t;
const dim_t DIM_MAX = std::numeric_limits<dim_t>::max();
const dim_t DIM_MIN = std::numeric_limits<dim_t>::min();

struct Point {
	Point() {};
	Point(dim_t x, dim_t y): x(x), y(y) {}
	dim_t x = 0;
	dim_t y = 0;
	bool marked = false;

	bool operator==(const Point& other) const {
		return this->x == other.x && this->y == other.y;
	}

	bool operator!=(const Point& other) const {
		return !this->operator==(other);
	}
};

struct Rect {
	Rect(dim_t x, dim_t y, dim_t width, dim_t height): x(x), y(y), width(width), height(height) {}
	dim_t x = 0;
	dim_t y = 0;
	dim_t width = 0;
	dim_t height = 0;
};

class Polygon : public vector<Point> {
public:
	Polygon() {}

	Polygon(std::initializer_list<Point> list) : vector<Point>(list) {
	}
	dim_t offsetx = 0;
	dim_t offsety = 0;
	dim_t x = 0;
	dim_t y = 0;
	dim_t width = 0;
	dim_t height = 0;
};

const Point INVALID_POINT(DIM_MAX,DIM_MAX);
const Rect INVALID_RECT(DIM_MAX,DIM_MAX,DIM_MAX,DIM_MAX);

const dim_t FLOAT_TOL = pow(10, -9);


bool _almostEqual(dim_t a, dim_t b, dim_t tolerance = FLOAT_TOL);
bool _withinDistance(const Point& p1, const Point& p2, dim_t distance);
dim_t _degreesToRadians(dim_t angle);
dim_t _radiansToDegrees(dim_t angle);
Point _normalizeVector(const Point& v);
bool _onSegment(const Point& A,const Point& B, const Point& p);
Point _lineIntersect(const Point& A,const Point& B,const Point& E,const Point& F, bool infinite = false);

namespace QuadradicBezier {

struct Segment {
	Point p1;
	Point p2;
	Point c1;
};

bool isFlat(const Point& p1, const Point& p2, const Point& c1, dim_t tol);
std::pair<Segment, Segment> subdivide(const Point& p1, const Point& p2,	const Point& c1, dim_t t);
list<Point> linearize(const Point& p1, const Point& p2, const Point& c1, dim_t tol);

} //QuadraticBezier

namespace CubicBezier {
struct Segment {
	Point p1;
	Point p2;
	Point c1;
	Point c2;
};

bool isFlat(const Point& p1, const Point& p2, const Point& c1, const Point& c2,	dim_t tol);
std::pair<Segment,Segment> subdivide(const Point& p1, const Point& p2, const Point& c1, const Point& c2, const dim_t& t);
list<Point> linearize(const Point& p1, const Point& p2, const Point& c1, const Point& c2, dim_t tol);

} //CubeBezier

namespace Arc {
struct CenterArc {
	Point center;
	dim_t rx;
	dim_t ry;
	dim_t theta;
	dim_t extent;
	dim_t angle;
};

struct SvgArc {
	Point p1;
	Point p2;
	dim_t rx;
	dim_t ry;
	dim_t angle;
	int8_t largearc;
	int8_t sweep;
};

SvgArc centerToSvg(const Point& center, const dim_t& rx, const dim_t& ry, dim_t theta1, const dim_t& extent, const dim_t& angleDegrees);
CenterArc svgToCenter(const Point& p1, const Point& p2, dim_t rx,	dim_t ry, const dim_t& angleDegrees, const int8_t& largearc, const int8_t& sweep);
list<Point> linearize(const Point& p1, const Point& p2, dim_t rx, dim_t ry,		dim_t angle, int8_t largearc, int8_t sweep, dim_t tol);

} //Arc

enum PointInPolygonResult {
	INSIDE,
	OUTSIDE,
	INVALID
};

Rect getPolygonBounds(const Polygon& polygon);
PointInPolygonResult pointInPolygon(const Point& point, const Polygon& polygon);
dim_t polygonArea(const Polygon& polygon);
bool intersect(const Polygon& A, const Polygon& B);
std::vector<Point> polygonEdge(const Polygon& polygon, Point normal);
dim_t pointLineDistance(const Point& p, const Point& s1, const Point& s2, Point normal, bool s1inclusive, bool s2inclusive);
dim_t pointDistance(const Point& p, const Point& s1, const Point& s2,		Point normal, bool infinite = false);
dim_t segmentDistance(const Point& A, const Point& B, const Point& E,	const Point& F, const Point& direction);
dim_t polygonSlideDistance(Polygon A, Polygon B, const Point& direction, bool ignoreNegative);
dim_t polygonProjectionDistance(Polygon A, Polygon B, const Point& direction);
// returns true if point already exists in the given nfp
bool inNfp(const Point& p, const std::vector<Polygon>& nfp);
Point searchStartPoint(Polygon A, Polygon B, bool inside,	const vector<Polygon>& NFP = {});
bool isRectangle(const Polygon& poly, dim_t tolerance = FLOAT_TOL) ;
// returns an interior NFP for the special case where A is a rectangle
vector<Polygon> noFitPolygonRectangle(const Polygon& A, const Polygon& B);

// given a static polygon A and a movable polygon B, compute a no fit polygon by orbiting B about A
	// if the inside flag is set, B is orbited inside of A rather than outside
	// if the searchEdges flag is set, all edges of A are explored for NFPs - multiple
vector<Polygon>	noFitPolygon(Polygon A, Polygon B, bool inside, bool searchEdges);

// given two polygons that touch at at least one point, but do not intersect. Return the outer perimeter of both polygons as a single continuous polygon
// A and B must have the same winding direction
Polygon polygonHull(Polygon A, Polygon B);
Polygon rotatePolygon(const Polygon& polygon, dim_t angle);
} //GeometryUtil
