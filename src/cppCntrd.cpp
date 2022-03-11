#include <iostream>
#include <vector>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <sstream>

extern "C"{
#include <fitsio.h>
}

/*
(420.08584608204205, 347.36249560529654)
*/

double sum2Dvec(const std::vector< std::vector<float>> &vec ){
    
    double sumVal = 0.0;
    int subFrameSize = vec[0].size();
    
    for (int ii = 0; ii < subFrameSize; ii++){
        for (int jj = 0; jj < subFrameSize; jj++){
            sumVal = vec[ii][jj] + sumVal;
        }
    }

    return sumVal;
}

std::vector< std::vector<float>> elementWise2DMult(const std::vector< std::vector<float>> &vec1,
                                                   const std::vector< std::vector<float>> &vec2){
    
    int subFrameSize = vec1[0].size();
    std::vector<std::vector<float>> out;
    out.resize(subFrameSize, std::vector<float>(subFrameSize, 0));

    for (int ii = 0; ii < subFrameSize; ii++){
        for (int jj = 0; jj < subFrameSize; jj++){
            out[ii][jj] = vec1[ii][jj] * vec2[ii][jj];
        }
    }

    return out;
}

void cntrd(const std::vector< std::vector<float>> &imageArray,
           const int &subFrameHalfSize,
           double * xc, double * yc){
    
    const int subFrameSize = 2 * (subFrameHalfSize);

    std::vector< std::vector<float>> yyii;
    yyii.resize(subFrameSize, std::vector<float>(subFrameSize, 0));
    
    std::vector< std::vector<float>> xxii;
    xxii.resize(subFrameSize, std::vector<float>(subFrameSize, 0));

    // # Step 2: Calculate the signal-weighted x and y moments
    // # (yi,xi) = indices((2*subR,2*subR))

    for(int ii = 0; ii < subFrameSize; ii++){
        for (int val=0; val < subFrameSize; val++){
            xxii[ii][val] = val;
            yyii[val][ii] = val;
        }
    }

    // yc = sum(yi*sf)/sum(sf)
    // xc = sum(xi*sf)/sum(sf)

    double sumSF = sum2Dvec(imageArray);
    std::vector<std::vector<float>> yisf = elementWise2DMult(yyii, imageArray);
    std::vector<std::vector<float>> xisf = elementWise2DMult(xxii, imageArray);

    *yc = sum2Dvec(yisf)/sumSF;
    *xc = sum2Dvec(xisf)/sumSF;
}

int main(int argc, char **argv){

    int status = 0;
    fitsfile *fptr;
    int naxis, bitpix;
    long naxisn[3];
    int nhdus;
    int subframeWidth = atoi(argv[2]);

    //==============================================//
    // load fits data
    //==============================================//
    fits_open_file(&fptr, argv[1], READONLY, &status);
    if(status) fits_report_error(stderr, status);
    fits_get_num_hdus(fptr, &nhdus, &status);
    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, naxis, naxisn, &status);
    fits_get_img_type(fptr, &bitpix, &status);

    std::vector<float> image(naxisn[0] * naxisn[1]);
    std::vector<long> fpixel(naxis, 1);

    std::ofstream outfile;
    std::stringstream fileLine;
    fileLine << std::setprecision(13) << std::fixed;
    outfile.open("cppCntrd.dat", std::ios::out | std::ios::trunc );

    //==============================================//
    // loop over all images in 3D fits file
    //==============================================//
    for(int ii = 0; ii < naxisn[2]; ii++){
        
        fpixel[2] = ii+1;
        fits_read_pix(fptr, TFLOAT, &fpixel[0], naxisn[0]*naxisn[1], NULL, &image[0], NULL, &status);
        
        //==============================================//
        // vector to 2d vector
        //==============================================//
        std::vector< std::vector<float>> imageArray;
        imageArray.resize(subframeWidth, std::vector<float>(subframeWidth, 0));
        for(unsigned int ii=0; ii < image.size(); ii++){
            imageArray[ii/subframeWidth][ii % subframeWidth] = image[ii];
        }

        //==============================================//
        // Try cntrd code
        //==============================================//
        double xc = 0;
        double yc = 0;
        cntrd(imageArray, (int) subframeWidth/2, &xc, &yc);

        // std::cout << "XC: " << xc << "\tYC: " << yc << "\n";
        fileLine << xc << ", " << yc << "\n";
    }
    outfile << fileLine.str();
    
    // std::cout << "Done" << std::endl;
    outfile.close();
    fits_close_file(fptr, &status);
    return 0;
}
