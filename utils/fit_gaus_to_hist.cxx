
#include <APEX_utils.h>
// ROOT headers
#include <TH1.h> 
#include <TF1.h>
#include <TAxis.h>  
#include <TFitResultPtr.h> 
#include <TFitResult.h>
#include <TString.h>
// stdlib headers
#include <stdexcept> 
#include <cmath> 

namespace APEX
{
namespace utils
{
    
/// @brief Attempts to fit a gaussian to a histogram, with a constant background
/// @param hist histogram to fit
/// @param radius radius of fit
/// @param center center of gaus, to be optimized by fit (value of arg. overwritten)
/// @param sigma sigma of gaus, to be optimized by fit (value of arg. overwritten) 
/// @param do_draw if true, then the fit will be drawn 
void fit_gaus_to_hist( 
    TH1 *hist, 
    double radius, 
    double &center, 
    double &sigma, 
    bool do_draw
) 
{
    //first guess for the center of the fit. we assume here that the gaus we wanna fit is the tallest peak in the hist. 
    center = hist->GetBinCenter( hist->GetMaximumBin() ); 
    
    //first guess for sigma 
    sigma = radius/4; 

    auto x_ax = hist->GetXaxis(); 
    
    double x_min = std::max( center-radius, x_ax->GetXmin() );    
    double x_max = std::min( center+radius, x_ax->GetXmax() );

    //the guess for the 'base' will be based on the bins on either extreme end of the fit
    double base = 0.5*(
        hist->GetBinContent( x_ax->GetBinCenter(x_min) ) +
        hist->GetBinContent( x_ax->GetBinCenter(x_max) )
    ); 

    //amplitude over background
    double amplitude = hist->GetMaximum() - base; 
    
    //our function to fit the histogram with   
    auto gaus_fit = new TF1( "gaus_fit", "gaus(0) + [3]", x_min, x_max ); 
  
    //set the parameters 
    gaus_fit->SetParameter( 0, amplitude ); 
    gaus_fit->SetParameter( 1, center ); 
    gaus_fit->SetParameter( 2, sigma ); 
    gaus_fit->SetParameter( 3, base ); 
    
    auto fit_ptr = hist->Fit("gaus_fit", (do_draw ? "L R S" : "L N R S Q")); 
    
    if (fit_ptr->IsValid()==false) {
        throw std::logic_error("<fit_gaus_to_hist>: fit failed."); 
        return; 
    }

    center = fit_ptr->Parameter(1); 
    sigma  = std::fabs(fit_ptr->Parameter(2)); 
    
    return; 
}

}
}