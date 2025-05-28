#include <cmath>
#include <vector>
#include <chrono>
using namespace std::chrono;
#include <iostream>
// #include <flann/flann.hpp>
#include <stdio.h>
#include "lbfgs.h"
#include <sstream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "std_image_write.h"

#define VOL_FLUID 0.6

#include <random>
static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform(0, 1);

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& a) {
    return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
Vector operator-=(Vector a, const Vector& b) {
    a[0] -= b[0];
    a[1] -= b[1];
    a[2] -= b[2];
    return a;
}
Vector operator+=(Vector a, const Vector& b) {
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
    return a;
}

class Polygon {  
public:
    std::vector<Vector> vertices;

    double area() {
        if (vertices.size() < 3) {return 0;}
        // return forumla from slides
        double s = 0.;
        for (int i=0; i<vertices.size(); i++) {
            int ip = (i==vertices.size()-1) ? 0 : (i + 1);
            s += vertices[i][0] * vertices[ip][1] - vertices[ip][0] * vertices[i][1];
        }
        return std::abs(s/2.);
    }

    Vector centroid() {
        if (vertices.size() < 3) {return Vector(0, 0, 0);}
        // return forumla from slides
        Vector c = Vector(0, 0, 0);
        for (int i=0; i<vertices.size(); i++) {
            int ip = (i==vertices.size()-1) ? 0 : (i + 1);
            double crossP = (vertices[i][0] * vertices[ip][1] - vertices[ip][0] * vertices[i][1]);
            c -= (vertices[i] + vertices[ip])*crossP; // because of minus sign, the orientation of the circles is wrong
        }
        double a = area();
        c = c/(6. * a);
        return c;
    }

    double integral_square_distance(const Vector& Pi) {
        if (vertices.size() < 3) {return 0;}
        double s = 0;
        // decompose polygon in n-2 triangles (n being the number of vertices)
        for (int t=1; t<vertices.size()-1; t++) {
            Vector c[3] = {vertices[0], vertices[t], vertices[t+1]}; //triangle

            double integralT = 0;
            for (int k=0; k<3; k++) {
                for (int l=k; l<3; l++) {
                    integralT += dot(c[k] - Pi, c[l] - Pi);
                }
            }
            Vector edge1 = c[1] - c[0];
            Vector edge2 = c[2] - c[0];
            // double areaT = 0.5 * std::abs((c[1][1]-c[0][1])*(c[2][0]-c[0][0]) - (c[1][0]-c[0][0])*(c[2][0]-c[0][0])); //formula cross product
            double areaT = 0.5 * std::abs((edge1[0]*edge2[1]) - (edge1[1]*edge2[0]));
            s += integralT * areaT / 6.;
        }
        return s;
    }
};  
     
// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon>& polygons, std::string filename, const std::vector<Vector>* points = NULL, std::string fillcol = "none") {
	FILE* f = fopen(filename.c_str(), "w+");
	fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
	for (int i = 0; i < polygons.size(); i++) {
		fprintf(f, "<g>\n");
		fprintf(f, "<polygon points = \"");
		for (int j = 0; j < polygons[i].vertices.size(); j++) {
			fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
		}
		fprintf(f, "\"\nfill = \"%s\" stroke = \"white\"/>\n", fillcol.c_str());
		fprintf(f, "</g>\n");
	}

	if (points) {
		fprintf(f, "<g>\n");		
		for (int i = 0; i < points->size(); i++) {
			fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r = \"3\" fill=\"white\" />\n", (*points)[i][0]*1000., 1000.-(*points)[i][1]*1000);
		}
		fprintf(f, "</g>\n");

	}

	fprintf(f, "</svg>\n");
	fclose(f);
}

int sgn(double x) {
    if (x > 0) return 0;
    if (x < 0) return 1;
    return 0;
}

void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W*H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W-1., std::max(0., W * bminx));
        bminy = std::min(H-1., std::max(0., H * bminy));
        bmaxx = std::max(W-1., std::max(0., W * bmaxx));
        bmaxy = std::max(H-1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
                    int sign = sgn(det);
                    if (prevSign == 0) prevSign = sign; else
                        if (sign == 0) sign = prevSign; else
                        if (sign != prevSign) {
                            isInside = false;
                            break;
                        }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                    double distEdge = std::abs(det)/ edgeLen;
                    double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                    if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {
                    //if (i < N) {   // the N first particles may represent fluid, displayed in blue
                    //  image[((H - y - 1)*W + x) * 3] = 0;
                    //  image[((H - y - 1)*W + x) * 3 + 1] = 0;
                    //  image[((H - y - 1)*W + x) * 3 + 2] = 255;
                    //}
                    if (mindistEdge <= 2) {
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 0;
                    }

                }
                
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

// For each point P_i
// construct Vor(P_i) # vornoi cell of the point P_I
// V = Square # a big square (initial cell)
// for each point P_j
        // V = cut(V, P_i, P_j)

/*
Polygon clipPolygon(Polygon subjectPolygon, Polygon clipPolygon) {
    Polygon outPolygon;
    for(const Vector& clipEdge : clipPolygon.vertices) { // For each edge of the clip polygon
        // Clip the subjectPolygon by a half-space
        for (int i = 0; i < subjectPolygon.vertices.size(); i++) { // For each vertex of the subject polygon
            // Test the subject polygon edge with vertices (i-1, i)
            Vector curVertex = subjectPolygon.vertices[i];
            Vector prevVertex = subjectPolygon.vertices[(i>0)?(i-1):(subjectPolygon.vertices.size()-1)];
            // Compute inter. between the infinite line supported by clipEdge and edge (i-1,i)
            Vector intersection = intersect(prevVertex, curVertex, clipEdge);
            if (curVertex inside clipEdge) {
                if (prevVertex not inside clipEdge) {
                    // The subject polygon edge crosses the clip edge, and we leave the clipping area
                    outPolygon.vertices.add(intersection);
                }
                outPolygon.vertices.add(curVertex);
            }
            else if (prevVertex inside clipEdge) {
                // The subject polygon edge crosses the clip edge, and we enter the clipping area
                outPolygon.vertices.add(intersection);
            }
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
}
*/
// TD2
// Semi discret optimal transprt (formula 1 and 2) + area of polygon + other

class VornoiDiagram {
public:
    // VornoiDiagram() {};

    VornoiDiagram() {
        N_disk = 100;
        unit_disk.resize(N_disk);
        for (int i = 0; i<N_disk; i++) {
            double theta = i*2. * M_PI / (double)N_disk;
            unit_disk[i] = Vector(cos(-theta), sin(-theta), 0);
        }
    };

    std::vector<Vector> unit_disk;
    int N_disk;

    Polygon clip_by_edge(const Polygon& V, const Vector &u, const Vector& v) {
        // Sutherman-Hodgman algo
        Vector N(v[1] - u[1], u[0] - v[0], 0); // check orientation otherwise, will keep all outside of polygon
        Polygon result;
        for (int i=0; i<V.vertices.size(); i++) {
            const Vector &A = V.vertices[(i == 0)? V.vertices.size()-1:i-1];
            const Vector &B = V.vertices[i];

            if (dot(u-B, N) >= 0) { // B inside = B closer to P0 than Pi
                if (dot(u-A, N) < 0) { // A outside
                    double t = dot(u-A, N) / dot(B-A, N);
                    Vector P = A + t*(B-A);
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else {
                if (dot(u-A, N) >= 0) { // A inside
                    // Vector M = (P0 + Pi) * 0.5;
                    double t = dot(u-A, N) / dot(B-A, N);
                    Vector P = A + t*(B-A);
                    result.vertices.push_back(P);
                }
            }
        }
        return result; // the clipped polygon
    }

    Polygon clip_by_bisector(const Polygon& V, const Vector &P0, const Vector& Pi, double w0, double wi) {
        Vector M = (P0 + Pi) * 0.5;
        Vector Mprime = M + ((w0 - wi) / (2*(P0-Pi).norm2()) * (Pi - P0)); // formula

        Polygon result;
        for (int i=0; i<V.vertices.size(); i++) {
            const Vector &A = V.vertices[(i == 0)? V.vertices.size()-1:i-1];
            const Vector &B = V.vertices[i];

            if ((B - P0).norm2() - w0 <= (B - Pi).norm2() - wi) { // B inside = B closer to P0 than Pi
                if ((A - P0).norm2() - w0 >= (A - Pi).norm2() - wi) { // A outside
                    double t = dot(Mprime-A,Pi-P0) / dot(B-A, Pi - P0);
                    Vector P = A + t*(B-A);
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else {
                if ((A - P0).norm2() - w0 <= (A - Pi).norm2() - wi) { // A inside
                    // Vector M = (P0 + Pi) * 0.5;
                    double t = dot(Mprime-A,Pi-P0) / dot(B-A, Pi - P0);
                    Vector P = A + t*(B-A);
                    result.vertices.push_back(P);
                }
            }
        }
        return result; // the clipped polygon
    }

    void compute() {
        Polygon square;
        square.vertices.push_back(Vector(0,0));
        square.vertices.push_back(Vector(0,1));
        square.vertices.push_back(Vector(1,1));
        square.vertices.push_back(Vector(1,0));

        cells.resize(points.size());

        /*
        for (auto& w : weights) {
            w = uniform(engine)* 0.1;
        }
        */

        #pragma omp parallel for schedule( dynamic, 1 )
        for (int i = 0; i < points.size(); i++) {
            Polygon V = square;
            for (int j=0; j< points.size(); j++) {
                if (i == j) continue;

                V = clip_by_bisector(V, points[i], points[j], weights[i], weights[j]);
            }
            // FLUID
            // clip V by disk
            for (int j=0; j<N_disk; j++) {
                double radius = sqrt(weights[i] - weights[weights.size()-1]);
                Vector u = unit_disk[j] * radius + points[i];
                Vector v = unit_disk[(j+1)%N_disk] * radius + points[i];
                V = clip_by_edge(V, u, v);
            }

            cells[i] = V;
        }
    }

    std::vector<Vector> points; // list of points
    std::vector<Polygon> cells;// list of polygons
    std::vector<double> weights; // = std::vector<double>(100);
};

class OptimalTransport {
public: 
    OptimalTransport() {};

    int ret ;

    void optimise();
    /*
    void optimise() {
        int N = vor.weights.size();
        lbfgsfloatval_t fx;
        std::vector<double> weights(N, 0);

        lbfgs_parameter_t param;
        // Initialise the parameters for the L-BFGS optimisation.
        lbfgs_parameter_init(&param);
        // param.lineserach = LBFGS_LINESEARCH_BACKTRACKING;

        int ret = lbfgs(N, &weights[0], &fx, evaluate, progress, (void*)this, &param);

        memcpy(&vor.weights[0], &weights[0], N * sizeof(weights[0]));
        vor.compute();
    }
    */

    VornoiDiagram vor;
};

/*
static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t* x,
    lbfgsfloatval_t* g,
    const int n,
    const lbfgsfloatval_t step
    )
{
    OptimalTransport* ot = (OptimalTransport*)(instance);

    memcpy(&ot->vor.weights[0], x, n * sizeof(x[0]));
    ot->vor.compute();

    int i;
    lbfgsfloatval_t fx = 0.0;

    for (i = 0;i < n;i += 1) {
        double current_area = ot->vor.cells[i].area();
        g[i] = current_area - 1./n;//- (formula from slides);

        fx += ot->vor.cells[i].integral_square_distance(ot->vor.points[i]) - x[i]*(current_area - 1./n);
    }
    return -fx;
}
*/

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t* x,
    lbfgsfloatval_t* g,
    const int n,
    const lbfgsfloatval_t step
)
{
    OptimalTransport* ot = (OptimalTransport*)(instance);

    memcpy(&ot->vor.weights[0], x, n * sizeof(x[0]));
    ot->vor.compute();

    lbfgsfloatval_t fx = 0.0;
    double sum_fluid_areas = 0.0;
    for (int i = 0; i < n; i++) {
        Polygon& cell = ot->vor.cells[i];
        if (cell.vertices.size() < 3) {
            g[i] = 0.0; // skip bad cells
            continue;
        }

        double area = cell.area();
        sum_fluid_areas += area;
        double integral = cell.integral_square_distance(ot->vor.points[i]);

        // g[i] = -(1.0 / n - area);  // == -∂g/∂w_i
        // Fluid
        g[i] = -(VOL_FLUID/n - area);
        fx += integral - x[i] * (area - 1.0 / n); // g(w)
        // FLUID
        fx += integral - x[i] * (area - VOL_FLUID / n); // g(w)
    }
    double estimated_air_volume = 1.0 - sum_fluid_areas;
    double desired_air_volume = 1.0 - VOL_FLUID;
    g[n-1] = -((1.0 - VOL_FLUID) - estimated_air_volume);
    fx += x[n-1] * (desired_air_volume - estimated_air_volume);

    return -fx; // because L-BFGS minimizes, but we maximize g
}


static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
    printf("Iteration %d:\n", k);
    // printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  fx = %f\n", fx);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

/*
void OptimalTransport::optimise() {
    int N = vor.weights.size();
    lbfgsfloatval_t fx;
    std::vector<double> weights(N, 0);
    // std::vector<double>& weights = vor.weights;

    lbfgs_parameter_t param;
    // Initialise the parameters for the L-BFGS optimisation.
    lbfgs_parameter_init(&param);
    // param.lineserach = LBFGS_LINESEARCH_BACKTRACKING;

    int ret = lbfgs(N, &weights[0], &fx, evaluate, progress, (void*)this, &param);

    memcpy(&vor.weights[0], &weights[0], N * sizeof(weights[0]));
    vor.compute();
}
*/

void OptimalTransport::optimise() {
    int N = vor.points.size(); 
    // vor.weights.resize(N); // match number of points
    std::vector<double> weights( N, 0 );
    memcpy(&weights[0], &vor.weights[0], N * sizeof(weights[0]));

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);

    lbfgsfloatval_t fx;
    int ret = lbfgs(N, &weights[0], &fx, evaluate, progress, (void*)this, &param);

    // vor.compute(); // final diagram with optimized weights
}

class Fluid {
public:
    Fluid(int N = 1000):N(N) {
        particles.resize(N);
        velocities.resize(N, Vector(0, 0, 0));
        for (int i=0; i<N; i++) {
            particles[i] = Vector(uniform(engine), uniform(engine), 0);
        }
        // fluid_volume = 0.6;

        ot.vor.points = particles;
        ot.vor.weights.resize(N+1);
        std::fill(ot.vor.weights.begin(), ot.vor.weights.end(), 1.);
        ot.vor.weights[N] = 0.99;
    } 

    void time_step(double dt) {

        double epsilon2 = 0.004 * 0.004;
        Vector g(0, -9.81, 0); //gravity
        double m_i = 100; //mass
        ot.vor.points = particles;
        ot.optimise();

        for (int i=0; i<particles.size(); i++) {
            Vector center_cell = ot.vor.cells[i].centroid();
            Vector spring_force = (center_cell - particles[i])/epsilon2;
            Vector all_forces = m_i*g + spring_force;
            velocities[i] += dt/m_i * all_forces;
            particles[i] += dt * velocities[i];

        }
    }

    void run_simulation() {
        double dt = 0.002;
        for (int i = 0; i< 100; i++) {
            time_step(dt);
            save_frame(ot.vor.cells, "test");
        }
    }

    OptimalTransport ot;
    std::vector<Vector> particles;
    std::vector<Vector> velocities;
    int N;
    double fluid_volume;
};

int main() {
    auto a = high_resolution_clock::now();

    int N = 1000;

    Fluid fluid = Fluid(100); // Fluid.fluid(100);
    fluid.run_simulation();
    exit(0);

    VornoiDiagram Vor;

    for (int i = 0; i < N; i++) {
        Vor.points.push_back(Vector(uniform(engine), uniform(engine), 0.0));
        Vor.weights.push_back( 0 );
    }
    Vor.compute();

    OptimalTransport ot;
    ot.vor = Vor;
    ot.optimise();

    save_svg(ot.vor.cells, "result.svg", &Vor.points); // , &Vor.points

    // End the timer
    auto b = high_resolution_clock::now();
    std::cout << "Took " << duration_cast<milliseconds>(b - a).count() << " milliseconds" <<  std::endl;
    return 0;
}