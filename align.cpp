#include "api.h"
#include "align.h"
#include "WorkLib.h"
#include <string>
#include <gtk/gtk.h>

using std::string;
using std::cout;
using std::endl;
using std::tie;
using std::make_tuple;
using std::tuple;


Subject::Subject(int a) : views(), type_change(a) {}

Subject::~Subject() {
}

void Subject::set_type(int t) {
    type_change = t;
}

void Subject::attach(Observer *o) {
    views.push_back(o);
}
void Subject::detach(Observer *o) {
    views.remove(o);
}
void Subject::notify() {
    std::list<Observer *>::iterator it;
    for (it = views.begin(); it != views.end(); it++) {
        (*it) -> update(type_change);
    }
}
//------------------------------------------

ImgModel::ImgModel() : Subject(0), image() {}

Image ImgModel::get_image() {
    return image;
}
void ImgModel::load(const char *p) {
    Image img = load_image(p);
    image = img;
    set_type(1);
    notify();
}
void ImgModel::make_filter(const vector<ImagePlugin *> &plugList, int numb_filt) {
    image = plugList[numb_filt] -> filter(image); 
    set_type(4);
    notify();
}

//------------------------------------------

Controller1::Controller1(ImgModel *i, ViewCons *c) : img(i), cons(c), manager() {}

Controller1::Controller1(const Controller1 &c) : img(c.img), cons(c.cons), manager(c.manager) {}
Controller1 &Controller1::operator=(const Controller1 &c) {
    img = c.img;
    cons = c.cons;
    manager = c.manager;
    return *this;
}

void Controller1::new_model(const char *p) {
    img -> load(p);
}
void Controller1::align_model() {
    img -> align(img -> get_image());
}
void Controller1::do_filter() {
    char pluginsDir[] = "plugins";
    vector<string> files;
    findLib(pluginsDir, files);
    vector<void *> pLib;
    loadLib(files, pLib);
    vector<string> nameFilt;
    loadPlugins(pLib, files, nameFilt, manager);
    int numb_filt = cons -> choose_filter(nameFilt);
    if (numb_filt >= 0) {
        img -> make_filter(manager.getPlugins(), numb_filt);
    }
}
//------------------------------------------

ViewCons::ViewCons(const ViewCons &view) : img(view.img), path(view.path) {}
ViewCons &ViewCons::operator= (const ViewCons &view) {
    img = view.img;
    path = view.path;
    return *this;
}
const char *ViewCons::get_path() {
    return path;
}


ViewCons::ViewCons(ImgModel *i, const char *p) : img(i), path(p) {
    img -> attach(this); 
}

ViewCons::~ViewCons () {
    img -> detach(this);
}

void ViewCons::update(int type) {

    if (type == 1) {
        write_info(type);
    } else if (type == 2) {
        write_info(type);
    } else if (type == 3) {
        write_info(type);
        save();
    } else if (type == 4) {
        write_info(type);
        save();
    }
}

void ViewCons::write_info(int type) {
    switch(type) {
        case 1 : {std::cout << "-- Image is load --" << endl; break;}
        case 2 : {std::cout << "-- Image is divided into color channels --" << endl; break;}
        case 3 : {std::cout << "-- Combination of channels is done --" << endl; break;}
        case 4 : {std::cout << "-- Filter applied --" << endl; break;}
        default: {}
    }
}

void ViewCons::save() {
    Image image = img -> get_image();
    save_image(image, path);
}

int ViewCons::choose_filter(const vector<string> &name_filters) {
    if (name_filters.size() > 0)
    {
        cout << endl << "Found filters:" << endl;
        for (size_t i=0; i<name_filters.size(); i++)
            cout << '[' << i << ']' << ' ' <<  (name_filters[i]) << endl;
        
        int pluginIndex = -1;
        while (pluginIndex >= signed (name_filters.size()) || pluginIndex < 0)
        {
            cout << endl << "Choose filter: ";
            cin >> pluginIndex;
            if (pluginIndex >= signed (name_filters.size()) || pluginIndex < 0)
            {
                cerr << "Invalid input: out of range " << "0-" << name_filters.size() - 1 << endl;
            }
        }
        return pluginIndex;
    } else {
        cout << "No filters found" << endl;
        return -1;
    }
}
//-------------------------------------------

struct pair_sh {
    
    float mse_min_1;
    float mse_min_2;
    int x_min_1;
    int y_min_1;
    int x_min_2;
    int y_min_2;


    pair_sh(float x1, int x2, int x3, float y1, int y2, int y3) : mse_min_1(x1), mse_min_2(y1), x_min_1(x2), y_min_1(x3),
        x_min_2(y2), y_min_2(y3) {}; 
};

pair_sh shift(Image img1, Image img2, Image img3, int h, int w, int x1, int y1, int x2, int y2, int delta1, int delta2) {
    int i, j;
    double mse_1, mse_2;
    double mse_min_1;
    double mse_min_2;
    int x_min_1, y_min_1;
    int x_min_2, y_min_2;
    uint r, g, b; 

    mse_1 = 0;
    mse_2 = 0;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            if (((i - x1) >= 0) && ((i - x1) < h) && ((j - y1) >= 0) && ((j - y1) <  w)) {
                b = std::get<2>(img1(i - x1,j - y1));
                g = std::get<1>(img2(i,j));
                mse_1 = mse_1 + (g - b)*(g - b);
                
            }
            if (((i - x2) >= 0) && ((i - x2) < h) && ((j - y2) >= 0) && ((j - y2) <  w)) {
                g = std::get<1>(img2(i,j));
                r = std::get<0>(img3(i - x2,j - y2));
                mse_2 = mse_2 + (g - r)*(g - r);
            }

        }
    }
    mse_1 = mse_1 / (h - abs(x1)) / (w - abs(y1));
    mse_2 = mse_2 / (h - abs(x2)) / (w - abs(y2));
    mse_min_1 = mse_1;
    mse_min_2 = mse_2;
    x_min_1 = x1;
    x_min_2 = x2;
    y_min_1 = y1;
    y_min_2 = y2;
    //printf("%6f %d %d \n", mse_min_1, x_min_1, y_min_1);
    //printf("%6f %d %d \n", mse_min_2, x_min_2, y_min_2);
    int step1, step2;
    int x1_beg = x1, y1_beg = y1;
    int x2_beg = x2, y2_beg = y2;
    for (step1 = -delta1; step1 <= delta1; step1++) {
        for (step2 = -delta2; step2 <= delta2; step2++) {
            x1 = x1_beg + step1;
            y1 = y1_beg + step2;
            x2 = x2_beg + step1;
            y2 = y2_beg + step2;
            //printf("%d %d ", x1, y1);
            mse_1 = 0;
            mse_2 = 0;
            for (i = 0; i < h; i++) {
                for (j = 0; j < w; j++) {
                    if (((i - x1) >= 0) && ((i - x1) < h) && ((j - y1) >= 0) && ((j - y1) <  w)) {
                        b = std::get<2>(img1(i - x1,j - y1));
                        g = std::get<1>(img2(i,j));
                        mse_1 = mse_1 + (g - b)*(g - b); 
                    }
                    if (((i - x2) >= 0) && ((i - x2) < h) && ((j - y2) >= 0) && ((j - y2) <  w)) {
                        g = std::get<1>(img2(i,j));
                        r = std::get<0>(img3(i - x2,j - y2));
                        mse_2 = mse_2 + (g - r)*(g - r);
                    }
                }
            }
            mse_1 = mse_1 / (h - abs(x1)) / (w - abs(y1));
            mse_2 = mse_2 / (h - abs(x2)) / (w - abs(y2));
            
            if (mse_1 < mse_min_1) {
                mse_min_1 = mse_1;
                x_min_1 = x1;
                y_min_1 = y1;
            }
            if (mse_2 < mse_min_2) {
                mse_min_2 = mse_2;
                x_min_2 = x2;
                y_min_2 = y2;
            }
            //printf("%6f\n", mse_1);
        }
    }
    //printf("%6f %d %d \n", mse_min_1, x_min_1, y_min_1);
    //printf("%6f %d %d \n", mse_min_2, x_min_2, y_min_2);

    pair_sh pair1(mse_min_1, x_min_1, y_min_1, mse_min_2, x_min_2, y_min_2);
    return pair1;

}

Image zoom2(Image img) {
    
    int hz = img.n_rows / 2;
    int wz = img.n_cols / 2;
    //printf("*%d %d*\n", hz, wz);
    int i, j;
    int sumr, sumg, sumb;
    Image img_zoom(hz, wz);
    for (i = 0; i < hz; i++) {
        for (j = 0; j < wz; j++) {
            sumr = std::get<0>(img(i * 2, j * 2)) + std::get<0>(img(i * 2, j * 2 + 1)) + std::get<0>(img(i * 2 + 1, j * 2)) + std::get<0>(img(i * 2 + 1, j * 2 + 1));
            sumr /= 4;
            sumg = std::get<1>(img(i * 2, j * 2)) + std::get<1>(img(i * 2, j * 2 + 1)) + std::get<1>(img(i * 2 + 1, j * 2)) + std::get<1>(img(i * 2 + 1, j * 2 + 1));
            sumg /= 4;
            sumb = std::get<2>(img(i * 2, j * 2)) + std::get<2>(img(i * 2, j * 2 + 1)) + std::get<2>(img(i * 2 + 1, j * 2)) + std::get<2>(img(i * 2 + 1, j * 2 + 1));
            sumb /= 4;
            img_zoom(i,j) = std::make_tuple(sumr, sumg, sumb);
        }
    }
    return img_zoom;

}

void ImgModel::align(Image srcImage)
{   
    int a = srcImage.n_rows / 3;
    int b = srcImage.n_cols;
    Image i1 = srcImage.submatrix(0, 0, a, b);
    Image i2 = srcImage.submatrix(a, 0, a, b);
    Image i3 = srcImage.submatrix(2 * a, 0, a, b);
 
    int height = srcImage.n_rows / 3;
    int width = srcImage.n_cols;
    int h = 0.9 * height;
    int w = 0.9 * width;
    Image img1 = srcImage.submatrix(0.05 * height, 0.05 * width, h, w);
    Image img2 = srcImage.submatrix(1.05 * height, 0.05 * width, h, w);
    Image img3 = srcImage.submatrix(2.05 * height, 0.05 * width, h, w);
    
    int i, j;
    int x_min_1, y_min_1;
    int x_min_2, y_min_2;
    
    bool isInterp = true;
    if (isInterp) {

        int size1 = img1.n_rows;
        int size2 = img1.n_cols;
        Image img_z1 = img1;
        Image img_z2 = img2;
        Image img_z3 = img3;
        int flag = 0;
        std::vector<Image> vect1, vect2, vect3;
        i = 0;
        while ((size1 > 400) && (size2 > 400)) {
            vect1.push_back(img_z1 = zoom2(img_z1));
            vect2.push_back(img_z2 = zoom2(img_z2));
            vect3.push_back(img_z3 = zoom2(img_z3));
            i++;
            
            size1 /= 2;
            size2 /= 2;
            //printf("%d %d %d\n", size1, size2, i);
            flag = 1;
        }
        
        if (flag) {
            i--;
            pair_sh pair1 = shift(vect1[i], vect2[i], vect3[i], size1, size2, 0, 0, 0, 0, 15, 15);
            
            //printf("%d %d\n", size1, size2);
            //pair_sh pair1 = shift(img1, img2, img3, h, w, 0, 0, 0, 0, 15, 15);
            x_min_1 = pair1.x_min_1; y_min_1 = pair1.y_min_1;
            x_min_2 = pair1.x_min_2; y_min_2 = pair1.y_min_2;
            //printf("%d %d | %d %d\n", x_min_1, y_min_1, x_min_2, y_min_2);

            for (j = i - 1; j >= 0; j--) {
                size1 *= 2;
                size2 *= 2;
                pair1 = shift(vect1[j], vect2[j], vect3[j], size1, size2, x_min_1 * 2, y_min_1 * 2, x_min_2 * 2, y_min_2 * 2  , 2, 2);

                x_min_1 = pair1.x_min_1; y_min_1 = pair1.y_min_1;
                x_min_2 = pair1.x_min_2; y_min_2 = pair1.y_min_2;
                //printf("%d %d | %d %d\n", x_min_1, y_min_1, x_min_2, y_min_2);
            }

            pair1 = shift(img1, img2, img3, h, w, x_min_1 * 2, y_min_1 * 2, x_min_2 * 2, y_min_2 * 2 , 2, 2);
            x_min_1 = pair1.x_min_1; y_min_1 = pair1.y_min_1;
            x_min_2 = pair1.x_min_2; y_min_2 = pair1.y_min_2;
            //printf("1) %d %d | %d %d\n", x_min_1, y_min_1, x_min_2, y_min_2);
            
        } else {
            pair_sh pair2 = shift(img1, img2, img3, h, w, 0, 0, 0, 0, 15, 15);

            x_min_1 = pair2.x_min_1; y_min_1 = pair2.y_min_1;
            x_min_2 = pair2.x_min_2; y_min_2 = pair2.y_min_2;

            //printf("2) %d %d | %d %d\n", x_min_1, y_min_1, x_min_2, y_min_2);
        }
    } else {
        pair_sh pair2 = shift(img1, img2, img3, h, w, 0, 0, 0, 0, 15, 15);
        x_min_1 = pair2.x_min_1; y_min_1 = pair2.y_min_1;
        x_min_2 = pair2.x_min_2; y_min_2 = pair2.y_min_2;
    }
//-------------------------------
    set_type(2);
    notify();
//-------------------------------
    
    int pr, pg, pb;
    Image img(a, b);

    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++) {
            pg = std::get<1>(i2(i, j));
            if ((i - x_min_1 >= 0) && (i - x_min_1 < a) && (j - y_min_1 >= 0) && (j - y_min_1 < b)) {
                pb = std::get<2>(i1(i - x_min_1, j - y_min_1));
            }
            if ((i - x_min_2 >= 0) && (i - x_min_2 < a) && (j - y_min_2 >= 0) && (j - y_min_2 < b)) {
                pr = std::get<0>(i3(i - x_min_2, j - y_min_2));
            }
            img(i,j) = std::make_tuple(pr, pg, pb);
        }
    }

//--------------------------------
    image = img;
    set_type(3);
    notify();
//--------------------------------

}
 //-----------------------------------------------------------------
/*class Unsharp {
private: Matrix<double> kernel = {{-1.0/6, -2.0/3, -1.0/6}, 
                               {-2.0/3, 13.0/3, -2.0/3}, 
                               {-1.0/6, -2.0/3, -1.0/6}};
public:
    tuple<uint, uint, uint> operator () (const Image &m) const
    {
        uint size = 2 * radius + 1;
        int r, g, b, new_r = 0, new_g = 0, new_b = 0;
        double r1 = 0, g1 = 0, b1 = 0;
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                tie(r, g, b) = m(i, j);                
                r1 += kernel(j, i) * r;
                g1 += kernel(j, i) * g;
                b1 += kernel(j, i) * b;
            }
        }
        new_r = r1;
        if (new_r > 255) new_r = 255;
        if (new_r < 0) new_r = 0;
        new_g = g1;
        if (new_g > 255) new_g = 255;
        if (new_g < 0) new_g = 0;
        new_b = b1;
        if (new_b > 255) new_b = 255;
        if (new_b < 0) new_b = 0;
        return make_tuple(new_r, new_g, new_b);
    }
    static const int radius = 1;

};


Image unsharp(Image src_image) {
    src_image = src_image.unary_map(Unsharp());
    return src_image;
}*/
//----------------------------------------------------------------------
/*class Gray_world {
private:
    float coef_r;
    float coef_g;
    float coef_b;
public:
    Gray_world(float r, float g, float b) : coef_r(r), coef_g(g), coef_b(b) {};

    tuple<uint, uint,uint> operator () (const Image &m) const { 
        uint r, g, b;
        tie(r, g, b) = m(0, 0);
        //printf("%4f %4f %4f)=> ", coef_r, coef_g, coef_b);
        r = r * coef_r;
        g = g * coef_g;
        b = b * coef_b;
        if (r > 255) {r = 255;};
        if (g > 255) {g = 255;};
        if (b > 255) {b = 255;};
        //printf("%d %d %d \n", r, g, b);
        return make_tuple(r, g, b);
    }
    static const int radius = 0;
};

Image gray_world(Image src_image) {

    uint h = src_image.n_rows;
    uint w = src_image.n_cols;
    long int r, g, b, sum_r = 0, sum_g = 0, sum_b = 0;
    for (uint i = 0.05 * h; i < 0.9 * h; ++i) {
        for (uint j = 0.05; j < 0.9 * w; ++j) {
            tie(r, g, b) = src_image(i, j);
            sum_r += r;
            sum_g += g;
            sum_b += b;
        }
    }
    double sum = sum_r + sum_g + sum_b;
    sum /= 3;
    src_image = src_image.unary_map(Gray_world(float (sum / sum_r), float (sum / sum_g), float (sum / sum_b)));
    return src_image;
}
*/
//---------------------------------------------------------------------


Image white_balance(Image src_image) {

    std::vector<int> red(256, 0), green(256, 0), blue(256, 0);
    uint i, j;
    uint h = src_image.n_rows;
    uint w = src_image.n_cols;

    Image image_tmp = src_image.deep_copy();
    image_tmp = image_tmp.submatrix( 0.05 * h, 0.05 * w, 0.9 * h, 0.9 * w);
    h = 0.9 * h;
    w = 0.9 * w;
    
    int r, g, b;
    float ravg = 0, gavg = 0, bavg = 0;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            tie(r, g, b) = image_tmp(i, j);
            red[r]++;
            green[g]++;
            blue[b]++;
            ravg += r;
            gavg += g;
            bavg += b;
        }
    }
    ravg /= h * w;
    gavg /= h * w;
    bavg /= h * w;
    float yavg = 0.299 * ravg + 0.587 * gavg + 0.114 * bavg;
    int ri = 0;
    int gi = 0;
    int bi = 0;

    for (i = 1; i <= 255; i++) {
        if (red[i] != 0) {
            red[i] += red[ri];
            ri = i;
        }
        if (green[i] != 0) {
            green[i] += green[gi];
            gi = i;
        }
        if (blue[i] != 0) {
            blue[i] += blue[bi];
            bi = i;
        }
    }
    int rmin = 0, gmin = 0, bmin = 0;
    while (red[rmin] == 0) {
        rmin++;  
    }
    while (green[gmin] == 0) {
        gmin++;
    }

    while (blue[bmin] == 0) {
        bmin++;
    }
    float x;
    for (i = 0; i <= 255; i++) {
        if (red[i] != 0) {
            x = (1.0 / (h * w - red[rmin]))*(red[i] - red[rmin]) * (255.0);
            red[i] = x;
        }
        if (blue[i] != 0) {
            x = (1.0 / (h * w - blue[bmin]))*(blue[i] - blue[bmin]) * (255.0);
            blue[i] = x;
        }
        if (green[i] != 0) {
            x = (1.0 / (h * w - green[gmin]))*(green[i] - green[gmin]) * (255.0);
            green[i] = x;
        }
    }
    Image image_hist = image_tmp.deep_copy();
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            tie(r, g, b) = image_tmp(i, j);
            image_hist(i, j) = make_tuple(red[r], green[g], blue[b]);
        }
    }

    float y;
    float cr, cb;
    Matrix<tuple<float, float, float>> ycb_img(h, w); 
    float y_f, cr_f, cb_f;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            tie(r, g, b) = image_hist(i, j);
            y_f = 0.299 * r + 0.587 * g + 0.114 * b;
            cr_f = (r - y_f) * 0.713;
            cb_f = (b - y_f) * 0.564;
            y = y_f;
            cr = cr_f;
            cb = cb_f;
            ycb_img(i, j) = make_tuple(y, cr, cb);
           
        }
    }

    
    std::vector<float> white0;
    std::vector<float> white1;
    std::vector<float> white2;
    uint count = 0;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            tie(y, cr, cb) = ycb_img(i, j);
            if ((y >= 210) && (cr <= 3) && (cr >= -3) && (cb <= 3) && (cb >= -3)) {
                white0.push_back(y);
                white1.push_back(cr);
                white2.push_back(cb);
                count++;
            }
        }
    }
    if (count != 0) {
        float ybri = -1, crbri = -4, cbbri = -4;
        float yav = 0, crav = 0, cbav = 0;
        for (i = 0; i < count; i++) {
            if (white0[i] > ybri) {
                ybri = white0[i];
                crbri = white1[i];
                cbbri = white2[i];
            } else if (abs(white0[i] - ybri) < 0.00001) {
                if (white1[i] + white2[i] < crbri + cbbri) {
                    crbri = white1[i];
                    cbbri = white2[i];
                }
            }
            yav += white0[i];
            crav += white1[i];
            cbav += white2[i];
        }
        yav /= count;
        crav /= count;
        cbav /= count;

        float ravg_h = yav + 1.402 * crav; 
        float gavg_h = yav - 0.34414 * cbav - 0.71414 * crav;
        float bavg_h = yav + 1.772 * cbav;
        float crl = std::min(crav, crbri);
        float cru = std::max(crav, crbri);
        float cbl = std::min(cbav, cbbri);
        float cbu = std::max(cbav, cbbri);

        count = 0;
        float rw = 0, gw = 0, bw = 0; 
        for (i = 0; i < h; i++) {
            for (j = 0; j < w; j++) {
                tie(y, cr, cb) = ycb_img(i, j);
                if ((y <= ybri) && (y >= yav) && (cr <= cru) && (cr >= crl) && (cb <= cbu) && (cb >= cbl) ) {
                    tie(r, g, b) = image_hist(i, j);
                    rw += r;
                    gw += g;
                    bw += b;
                    count ++;
                }
            }
        }
        rw /= (count);
        gw /= (count);
        bw /= (count);
        float yw = 0.299 * rw + 0.587 * gw + 0.114 * bw;
        float rscale = float (yw / rw);
        float gscale = float (yw / gw);
        float bscale = float (yw / bw);
        float r_gwa = float (yavg / ravg);
        float g_gwa = float (yavg / gavg);
        float b_gwa = float (yavg / bavg);

        float rfact = 0, gfact = 0, bfact = 0;
        int flag = 1;

        if (((bavg_h + 3) >= gavg_h) && (bavg_h > ravg_h)) {
            rfact = rscale;
            gfact = gscale;
            bfact = b_gwa;
        }
        else if (gavg_h + 3 > bavg_h) {
            rfact = rscale;
            gfact = g_gwa;
            bfact = bscale;
        }
        else if ((ravg_h > gavg_h) && (gavg_h > bavg_h)) {
            rfact = r_gwa;
            gfact = gscale;
            bfact = bscale;
        } else { flag = 0; }
        if (flag) {
            for (i = 0; i < src_image.n_rows; i++) {
                for (j = 0; j < src_image.n_cols; j++) {
                    tie(r, g, b) = src_image(i, j);
                    r = float (rfact * r);
                    g = float (gfact * g);
                    b = float (bfact * b);
                    if (r > 255) {r = 255;};
                    if (g > 255) {g = 255;};
                    if (b > 255) {b = 255;};
                    src_image(i, j) = make_tuple(r, g, b);
                }
            }
        }
    } 

    return src_image;
}
Image custom(Image src_image, Matrix<double> kernel) {
    // Function custom is useful for making concrete linear filtrations
    // like gaussian or sobel. So, we assume that you implement custom
    // and then implement other filtrations using this function.
    // sobel_x and sobel_y are given as an example.
    return src_image;
}
//---------------------------------------------------------

Image autocontrast(Image src_image, double fraction) {
    int i, j;
    int h = src_image.n_rows;
    int w = src_image.n_cols;
    int r, g, b, y;
    std::vector<float> osv;
    float ymin = 255, ymax = 0;
    for (i = 0.05 * h; i < 0.9 * h; i++) {
        for (j = 0.05 * w; j < 0.9 * w; j++) {
            tie(r, g, b) = src_image(i, j);
            y = 0.2125 * r + 0.7154 * g + 0.0721 * b;
            osv.push_back(y);
        }
    }
    sort(osv.begin(), osv.end());
    ymin = osv[fraction * osv.size()];
    ymax = osv[osv.size() - fraction * osv.size() - 1];
    float d = ymax - ymin;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            tie(r, g, b) = src_image(i, j);
            r = (float (255 / d)) * (r - ymin);
            g = (float (255 / d)) * (g - ymin);
            b = (float (255 / d)) * (b - ymin);
            if (r < 0) { r = 0; };
            if (r > 255) { r = 255; };
            if (g < 0) { g = 0; };
            if (g > 255) { g = 255; };
            if (b < 0) { b = 0; };
            if (b > 255) { b = 255; };
            src_image(i, j) = make_tuple(r, g, b);
        }
    }
    return src_image;
}

//-------------------------------------------------------

Image mirror(Image src_image, int radius) {
    int h = src_image.n_rows + 2 * radius;
    int w =  src_image.n_cols + 2 * radius;
    int r = radius;
    int i, j;
    Image img(h, w);
    for (i = radius; i < h - radius; i++) {
        for (j = radius; j < w - radius; j++) {
            img(i, j) = src_image(i- r, j - r);
        }
    }
    for (i = r; i < h - r; i++) {
        for (j = 0; j < r; j++) {
            img(i, j) = src_image(i - r, j + 2 * (r - j) - 1 - r);
            img(i, w - r + j) = src_image(i - r, w - r + j - 2 * j - 1 - r);
        }
    }
    for (j = r; j < w - r; j++) {
        for (i = 0; i < r; i++) {
            img(i, j) = src_image(i + 2 * (r - i) - 1 - r, j - r);
            img(h - r + i, j) = src_image(h - r + i - 2 * i - 1 - r , j - r);
        }
    }
    for (j = 0; j < r; j++) {
        for (i = 0; i < r; i++) {
            img(i, j) = img(i + 2 * (r - i) - 1, j);
            img(h - r + i, j) = img(h - r + i - 2 * i - 1 , j);
        }
    }
    for (j = w - r; j < w; j++) {
        for (i = 0; i < r; i++) {
            img(i, j) = img(i + 2 * (r - i) - 1, j);
            img(h - r + i, j) = img(h - r + i - 2 * i - 1 , j);
        }
    }
    return img;
}
//----------------------------------------------------------------
class MedianF {
public:
    MedianF (int a) : radius(a) {};
    tuple<uint, uint, uint> operator () (const Image &m) const
    {
        uint size = 2 * radius + 1;
        uint r, g, b;
        std::vector<uint> r_pix, g_pix, b_pix;
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                tie(r, g, b) = m(i, j);
                r_pix.push_back(r);
                g_pix.push_back(g);
                b_pix.push_back(b);
            }
        }
        std::sort(r_pix.begin(), r_pix.end());
        std::sort(g_pix.begin(), g_pix.end());
        std::sort(b_pix.begin(), b_pix.end());
        r = r_pix[(size * size - 1) / 2];
        g = g_pix[(size * size - 1) / 2];
        b = b_pix[(size * size - 1) / 2];

        return make_tuple(r, g, b);
    }
    int radius;
};

Image median(Image src_image, int radius) {
    Image img = mirror(src_image, radius);
    img = img.unary_map(MedianF(radius));
    src_image = img.submatrix(radius, radius, src_image.n_rows, src_image.n_cols);
    return src_image;
}
//------------------------------------------------------------------
Image median_linear(Image src_image, int radius) {
    Image img = mirror(src_image, radius);
    int rh[256] = {0};
    int gh[256] = {0};
    int bh[256] = {0};
    int h = img.n_rows;
    int w = img.n_cols;
    int rad = radius;
    Image res(h, w);
    int i, j, k1, k2, r, g, b, sum_r, sum_g, sum_b;
    for (i = rad; i < h - rad; i++ ) {
        for (k1 = 0; k1 < 256; k1++) {
            rh[k1] = 0;
            gh[k1] = 0;
            bh[k1] = 0;
        }
        for (k1 = -rad; k1 <= rad; k1++) {
            for (k2 = 0; k2 <= 2 * rad; k2++) {
                tie(r, g, b) = img(i + k1, k2);
                rh[r]++;
                gh[g]++;
                bh[b]++;
            }
        }
        for (j = rad; j < w - rad; j++) {
            if (j != rad) {
                for (k1 = -rad; k1 <= rad; k1++) {
                    tie(r, g, b) = img(i + k1, j - rad - 1);
                    rh[r]--;
                    gh[g]--;
                    bh[b]--;
                    tie(r, g, b) = img(i + k1, j + rad);
                    rh[r]++;
                    gh[g]++;
                    bh[b]++;
                }
            }
            sum_r = 0;
            sum_g = 0;
            sum_b = 0;
            k1 = -1;
            while ( sum_r < (2 * rad + 1) * (2 * rad + 1) / 2) {
                k1++;
                sum_r += rh[k1];
            }
            r = k1;
            k1 = -1;
            while ( sum_g < (2 * rad + 1) * (2 * rad + 1) / 2) {
                k1++;
                sum_g += gh[k1];
            }
            g = k1;
            k1 = -1;
            while ( sum_b < (2 * rad + 1) * (2 * rad + 1) / 2) {
                k1++;
                sum_b += bh[k1];
            }
            b = k1;
            res(i, j) = make_tuple(r, g, b);
        }
    }
    res = res.submatrix(rad, rad, src_image.n_rows, src_image.n_cols);
    return res;
}
//--------------------------------------------------------------------

