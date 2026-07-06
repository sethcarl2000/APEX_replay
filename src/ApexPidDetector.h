#ifndef ApexPidDetector_H
#define ApexPidDetector_H

#include <utility>
#include <vector> 
#include <ROOT/RVec.hxx>
#include <ApexVDCTrack.h>
//#include <ApexPidManager.h> 
class ApexPidManager; 

class ApexPidDetector {
public:   

    struct Cell {
        int row,col, id; 
        
        bool operator==(const Cell& rhs) const { return rhs.id==id; }
        bool operator<(const Cell& rhs) const { 
            if (row < rhs.row) return true; 
            if (col < rhs.col) return true;
            return false;  
        }
    };

private: 
    // we declare the ApexPidManager as a friend, so it can access the following private data members 
    friend class ApexPidManager;     

    // generate cells in the ususal order 
    std::vector<Cell> GenerateCells(int n_rows, int n_cols); 

    // calorimeter cells
    std::vector<Cell> fCells;  

    double fZ; //< z-position of the plane 

    int fN_rows; //< n. rows

    int fN_cols; //< n. cols

    double fCell_x0; //< center of first cell 'x' (m)
    double fCell_y0; //< center of first cell 'y' (m)

    double fCell_width_x; //< width of cells in x direction (m)
    double fCell_width_y; //< width of cells in y direction (m)

    const Cell* GetCell(int row, int col) const; 

public: 

    ApexPidDetector() {}; 

    virtual ~ApexPidDetector() {}; 

    /// @brief given array of ADC data, return the data of cell indexed by 'row / col'.  
    /// @param row row of cell
    /// @param col column of cell
    /// @param data array of ADC data
    /// @return adc data in cell [row, col]. returns '0' if cell doesn't exist  
    double GetVal(int row, int col, const ROOT::RVec<double>& data) const; 

    /// @brief Get the closest cell associated with this track 
    /// @param track the track to project 
    /// @return the cell the track intercepts with 
    const Cell* GetNearestCell(const ApexVDC::Track& track) const; 

    inline double GetZ() const { return fZ; }

    void GetCellXY(int row, int col, double& x, double& y) const; 

    const Cell* GetNearestCell(double x, double y, double dxdz, double dydz) const;  
}; 


#endif