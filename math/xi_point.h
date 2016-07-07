#ifndef XI_POINT_H
#define XI_POINT_H

//class to store xi points, to help send data from the Mesh to the VisualizationWindow
class xiPoint
{
public:
    unsigned x, y;  //coordinates (discrete)
    int zero, one, two;  //multiplicities of xi_0, xi_1, and xi_2 at this point ---- TODO: maybe should be unsigned?

    xiPoint(unsigned xc, unsigned yc, int m0, int m1, int m2);
    template <class Archive>
    void cerealize(Archive &archive) {
        archive(x, y, zero, one, two);
    }
};
#endif // XI_POINT_H
