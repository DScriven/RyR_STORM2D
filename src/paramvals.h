#ifndef PARAMVALS_H
#define PARAMVALS_H
#include <QRect>

struct ParamVals
{
    QRect  paramGeometry;
    double MinTetramerArea;
    double NeighbourhoodLimit;
    double Alpha1stPass;
    double LogDensityThreshold;
    uint MinimumBlinksPerCluster;
    bool SaveClusterData;
    bool SaveClusterBlinks;
    bool DispNoThreshold;
    bool DispExcludedData;
    bool DispNumbers;
    bool DispConvexHull;
    bool DispMinEllipse;
    bool DispBoundary;
};

#endif // PARAMVALS_H
