#include <ApexPidDetector.h>
#include <algorithm>
#include <TError.h> 
#include <limits> 
#include <stdexcept> 
#include <sstream> 
#include <ApexUtils.h> 

//_____________________________________________________________________________________________________________________
const ApexPidDetector::Cell* ApexPidDetector::GetCell(int row, int col) const
{
    auto it = std::find_if(
        fCells.begin(),
        fCells.end(),
        [row,col](const Cell& rhs){ return (rhs.row==row) && (rhs.col==col); }
    );

    if (it != fCells.end()) { return &(*it); }
    return nullptr;     
}
//_____________________________________________________________________________________________________________________
double ApexPidDetector::GetVal(int row, int col, const ROOT::RVec<double>& data) const
{
    const Cell* cell = GetCell(row, col);
    //check if this cell was found
    if ((int)data.size() != fN_rows*fN_cols) {
        std::ostringstream oss; 
        oss << "in <ApexPidDetector::"<<__func__<<"> size of data input ("<<data.size()<<") does not match number of cells ("<<fN_cols*fN_rows<<")"; 
        throw std::logic_error(oss.str());  
        return APEX::kNaN_double; 
    }
    if (cell) { return data[cell->id]; } else { return 0.; }
}
//_____________________________________________________________________________________________________________________
const ApexPidDetector::Cell* ApexPidDetector::GetNearestCell(const ApexVDC::Track& track) const
{
    //project the track onto the z-plane
    return GetNearestCell(track.FP_x(), track.FP_y(), track.FP_dx_dz(), track.FP_dy_dz());
}
//_____________________________________________________________________________________________________________________
const ApexPidDetector::Cell* ApexPidDetector::GetNearestCell(double x0, double y0, double dxdz, double dydz) const
{
    //change from focal-plane to transport coordiantes
    dxdz = dxdz + x0/6.; 
    
    //project the track onto the z-plane.
    double x = x0 + (dxdz*fZ); 
    double y = y0 + (dydz*fZ);

    //get the closest cell 
    int col = (int)std::round( (fCell_x0 - x)/fCell_width_x ); 
    int row = (int)std::round( (fCell_y0 - y)/fCell_width_y );

    if (col >= fN_cols) { col=fN_cols-1; } else { if (col < 0) col=0; }
    if (row >= fN_rows) { row=fN_rows-1; } else { if (row < 0) row=0; }

    return GetCell(row, col); 
}
//_____________________________________________________________________________________________________________________
std::vector<ApexPidDetector::Cell> ApexPidDetector::GenerateCells(int n_rows, int n_cols)
{
    fCells.clear(); fCells.reserve(n_rows*n_cols); 

    int id=0; 
    for (int row=0; row<n_rows; row++) 
        for (int col=0; col<n_cols; col++) 
            fCells.push_back({ row, col, id++ }); 
        
    return fCells; 
}
//_____________________________________________________________________________________________________________________
void ApexPidDetector::GetCellXY(int row, int col, double& x, double& y) const
{
    x = fCell_x0 - ((double)col)*fCell_width_x; 
    y = fCell_y0 - ((double)row)*fCell_width_y; 
}
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________