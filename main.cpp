#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <list>

using std::string;
using std::stringstream;
using std::cout;
using std::cerr;
using std::endl;
using std::numeric_limits;

#include "api.h"
#include "align.h"

//-------------------------------------------
int main(int argc, char const *argv[])
{

    ImgModel *img = new ImgModel;   
    ViewCons *view1 = new ViewCons(img, argv[2]);

    std::cout << view1 -> get_path() << "\n";
    Controller1 *contr = new Controller1(img, view1);
    contr -> new_model(argv[1]);
    contr -> align_model();
    if (argc == 4) {
        std::string filter(argv[3]);
        if (filter == "--filter") {
            contr -> do_filter();
        }
    }
    return 0;
}