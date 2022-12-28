#ifndef TETRAMER_H
#define TETRAMER_H
struct Tetramer
{
    int ClusterNo;
    float xcentre;
    float ycentre;
    int angle;
    char classification;

    Tetramer(int k, float x, float y, int ar, char c) : ClusterNo(k), xcentre(x), ycentre(y), angle(ar), classification(c) {}
    Tetramer() {}

    bool operator < (const Tetramer& tr) const
    {
        return (ClusterNo < tr.ClusterNo);
    }
};

#endif // TETRAMER_H
