#include "ncorr.h"
#include "DIC_3D.h"
#define OUTPUT_GENERATE_VIDEO

using namespace ncorr;

int main(int argc, char *argv[]) {
	if (argc != 2) {
		throw std::invalid_argument("Must have 1 command line input of either 'calculate' or 'load'");	
	}

	DIC_3D_Input dic_3d_input("3D_DIC_input.txt");
	dic_3d_input.debug_print(std::cout);

	std::string path_prefix = "";
	// Initialize DIC and strain information ---------------//
	DIC_analysis_input DIC_input;
	DIC_analysis_output DIC_output;

	// Determine whether or not to perform calculations or 
	// load data (only load data if analysis has already 
	// been done and saved or else throw an exception).
	std::string input(argv[1]);
	if (input == "load") {
		// Load inputs
		DIC_input =    		DIC_analysis_input::load(path_prefix + "save/DIC_input.bin");
		DIC_output = 	   DIC_analysis_output::load(path_prefix +"save/DIC_output.bin");
	} else if (input == "calculate") {
		// Set images
		const size_t images_count = 2;
		std::vector<Image2D> imgs;
		for (int i = 1; i <= images_count; ++i) {
		    std::ostringstream ostr;
		    ostr << path_prefix + "images/00" << std::setfill('0') << std::setw(2) << i << ".jpg";
		    imgs.push_back(ostr.str());
		}
		
		// Set DIC_input
		DIC_input = DIC_analysis_input(imgs, 							// Images
				               ROI2D(Image2D(path_prefix + "images/roi_full.png").get_gs() > 0.5),		// ROI
					       4,                                         		// scalefactor
					       INTERP::CUBIC_KEYS_PRECOMPUTE,			// Interpolation
					       SUBREGION::SQUARE,					// Subregion shape
					       35,     // in pixels!                                   		// Subregion radius
					       8,                                         		// # of threads
					       DIC_analysis_config::NO_UPDATE,				// DIC configuration for reference image updates
					       false);//true);							// Debugging enabled/disabled

		// Perform DIC_analysis    
		DIC_output = DIC_analysis(DIC_input);

		// Convert DIC_output to Eulerian perspective
		// DIC_output = change_perspective(DIC_output, INTERP::QUINTIC_BSPLINE_PRECOMPUTE);

		// Set units of DIC_output (provide units/pixel)
		DIC_output = set_units(DIC_output, "mm", 1.0);
		
		// Save outputs as binary
                // save(DIC_input, path_prefix + "save/DIC_input.bin");
                // save(DIC_output, path_prefix + "save/DIC_output.bin");
	} else {
		throw std::invalid_argument("Input of " + input + " is not recognized. Must be either 'calculate' or 'load'");	
	}	

	std::cout<<"disps:"<<DIC_output.disps.size()<<std::endl;
        

#ifdef OUTPUT_GENERATE_VIDEO
        save_DIC_video(path_prefix + "video/test_v_eulerian.avi", 
                       DIC_input, 
                       DIC_output, 
                       DISP::V,
                       0.5,		// Alpha		
                       15);		// FPS

        save_DIC_video(path_prefix + "video/test_u_eulerian.avi", 
                       DIC_input, 
                       DIC_output, 
                       DISP::U, 
                       0.5,		// Alpha
                       15);		// FPS

#endif
		
						  
		save_ncorr_data_over_img(path_prefix + "images/test_u_last_eulerian.jpg",
			DIC_input.imgs.back(),
			DIC_output.disps.back().get_u(),
			0.5,
			-600,
			600, //grade
			true,
			true,
			true,
			DIC_output.units,
			DIC_output.units_per_pixel,
			50,
			1.0,
			11,
			cv::COLORMAP_JET);

		save_ncorr_data_over_img(path_prefix + "images/test_v_last_eulerian.jpg",
			DIC_input.imgs.back(),
			DIC_output.disps.back().get_v(),
			0.5,
			0.0,
			0.6, //grade
			true,
			true,
			true,
			DIC_output.units,
			DIC_output.units_per_pixel,
			50,
			1.0,
			11,
			cv::COLORMAP_JET);

  	return 0;
}
