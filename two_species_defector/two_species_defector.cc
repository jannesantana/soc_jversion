// 
//  To compile: c++ abs_anti_mips_spinoff_ll_v2.cc cokus3.c  -o antimips -lgsl -lgslcblas -lm
// If in mac : c++ -Wall -I/opt/homebrew/Cellar/gsl/2.7.1/include/ -L/opt/homebrew/Cellar/gsl/2.7.1/lib/ abs_anti_mips_spinoff_ll_v2.cc cokus3.c -o antimips -lgsl -lgslcblas -lm
// using namespace std;


#include <sstream>
#include <map>
#include <iostream>
// #include <string.h>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <ctime>
#include <random>
#include <algorithm> 
#include <vector>
#define PI 3.14159265358979323846

// Constants and definitions 
double box; // Size of the 2D box (assuming a square box)
double CUTOFF;    // Interaction cutoff distance
int N_particles;  // Number of particles
int Cells = box / CUTOFF; // Number of cells per dimension
int T ;
double Dt;
double Vo;
double eta;
double R;
double sigma;
double epsilon;
double alig_str;
double acceptance_rate=9;
double op_now;
std::vector<int> particles_w_high_order;


std::map<std::string, double> readParams(const std::string& filename) {
    std::ifstream file(filename);
    std::map<std::string, double> params;
    std::string line, param_name;
    double param_value;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            if (!(iss >> param_name >> param_value)) { break; } // Error in line format
            params[param_name] = param_value;
        }
        file.close();
    } else {
        std::cerr << "Unable to open parameter file: " << filename << std::endl;
        exit(1);
    }
    return params;
}

// Function to create the filename based on parameters
std::string createFileName(const std::map<std::string, double>& params, const std::string& baseName) {
    std::ostringstream oss;
    oss << baseName;
    for (const auto& param : params) {
        oss << "_" << param.first << "_" << param.second;
    }
    oss << ".dat";
    return oss.str();
}


// int N;
	
double TopeX;
double TopeY;
	
int TempsInitial;
int TempsTotal;
// double Dt; 
int t_XY_Save=1;	




double atu_dist;

double dist;





struct Particle {
    double x, y; // Position
    double theta; //angle 
    double fx, fy; // colisional force 
    double alignment; 
    int neighbors; 
    double avg_ang_region;
    double defector;
    int if_defector;
    int cellIndex; // Cell index in the linked list
    double sd_t;
    double prob_of_accepting;
    std::vector<int> neigh_loc;
    double ini_posx, ini_posy;
    double aux_posx, aux_posy;
    double op_in_region;
};

double applyPBC(double coord) {
    if (coord >= box) {coord = coord - box;}
    if (coord < 0) {coord = coord + box;}
     return coord;
}

 
// angular boundary conditions  
double applyPBCangle(double angle) {
    if (angle >= 2*PI) return angle - 2*PI;
    if (angle < 0) return angle + 2*PI;
    return angle;
}

// associate position to cell 
int getCellIndex(double x, double y) {

    int cellX =  static_cast<int>(x / CUTOFF) % Cells; 
    int cellY  = static_cast<int>(y / CUTOFF) % Cells;

    if (cellX < 0) cellX += Cells; // negative correction
    if (cellY < 0) cellY += Cells;
    if (cellX < 0 || cellX >= Cells || cellY < 0 || cellY >= Cells) { //  in case shit goes wrong 
        
        std::cerr << "Cell index out of bounds: (" << cellX << ", " << cellY << ")" <<  std::endl;
        exit(1);
    }
    return cellY * Cells + cellX;
    
}


   std::random_device rd{}; 
    std::mt19937 gen(rd()); 
    std::normal_distribution<double> dis{0.0, 1.0};// Standard mersenne_twister_engine seeded with rd()
 ; 
    std::uniform_real_distribution<double> dis2(0, 1.0);


void initializeParticles(Particle* particles) {

    
    
    for (int i = 0; i < N_particles; ++i) {

        particles[i].x = dis2(gen) * box; // random initial consitions 
        particles[i].y = dis2(gen) * box;

        particles[i].theta = dis(gen)*PI;
		particles[i].fx = 0.0;
        particles[i].fy = 0.0;
        particles[i].alignment = 0.0;
		 particles[i].ini_posx = 0.0;
         particles[i].neighbors = 0;
        particles[i].ini_posy = 0.0;
        particles[i].avg_ang_region=0.0;
        particles[i].defector=0.0;
        // particles[i].if_defector= -1;
        particles[i].aux_posx=0.0;
        particles[i].aux_posy=0.0;
        particles[i].op_in_region = 0.0;
        particles[i].cellIndex = getCellIndex(particles[i].x, particles[i].y); // get cell index 

	}

	}

// void savePositions(const Particle* particles) { 
//     std::ofstream file;
   
//     file.open("particle_positions.dat", std::ios::app); // open the file 

//     // file << "Timestep " << timestep << "\n";
//     for (int i = 0; i < N_particles; ++i) {
//         file << particles[i].x << " " << particles[i].y  << " " << particles[i].theta <<"\n";
//     }

    

//     file.close();
// }

void savePositions(const Particle* particles, const std::string& file_name) {
    std::ofstream file;
    file.open(file_name, std::ios::app);
    for (int i = 0; i < N_particles; ++i) {
        file << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
    }
    file.close();
}


void buildLinkedList(Particle* particles, int* head, int* linkedList) {
    // Initialize head and linked list
    for (int i = 0; i < Cells * Cells; ++i) {
        // std::cout << "Cell " << i << ": ";
        head[i] = -1;
    }
    for (int i = 0; i < N_particles; ++i) {
        linkedList[i] = -1;
    }

    for (int i = 0; i < N_particles; ++i) { // initialize the complete linked list 
        int cellIndex = particles[i].cellIndex;
        linkedList[i] = head[cellIndex];
        head[cellIndex] = i;
    }

    for (int i = 0; i < Cells * Cells; ++i) { // neighbouring linked list 
        // std::cout << "Cell " << i << ": ";
        int j = head[i];
        
    }

   

}

void simulateInteractions(Particle* particles, int* head,  int* linkedList) {

    particles_w_high_order.clear();
    particles_w_high_order.shrink_to_fit();
    for (int i = 0; i < N_particles; ++i) {
        
        particles[i].fx = 0.0; // Reset forces
        particles[i].fy = 0.0;
        particles[i].alignment= 0.0;
        particles[i].neighbors = 1;
        particles[i].avg_ang_region=0.0;
        particles[i].defector=0.0;
        particles[i].neigh_loc.clear();
        particles[i].neigh_loc.shrink_to_fit();
        particles[i].op_in_region = 0.0;
        // particles[i].if_defector= -1;
        // particles[i].force_pili_x = 0.0;
        // particles[i].force_pili_y = 0.0;
       
    }

double order_in_region_x = 0.0; 
double order_in_region_y = 0.0; 

// double op_in_region= 0.0;
// int neigh_in_region;

 for (int i = 0; i < N_particles; ++i) {
   
Particle& pi = particles[i];
        int cellX = static_cast<int>(pi.x / CUTOFF)% Cells;
        int cellY = static_cast<int>(pi.y / CUTOFF)% Cells;
        // std::cout <<"cell index for particle at the moment" << std::endl;

        // Loop over neighboring cells
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                int neighborX = (cellX + dx + Cells) % Cells;
                int neighborY = (cellY + dy + Cells) % Cells;
                int neighborCellIndex = neighborY * Cells + neighborX; // pick up neighboring cell index 
                // if (cellX==0 & cellY==0) {
                //     std::cout << "neighbouring cells = (" << neighborX <<","<< neighborY <<")" <<std::endl;
                // }

                //  std::cout <<" searching for neighbors" << std::endl;

                int j = head[neighborCellIndex];
                while (j != -1) {
                    
                    
                    if (i != j) {
                        Particle& pj = particles[j];
                        double dx = pi.x - pj.x;
                        double dy = pi.y - pj.y;

                       

                       

                        // PBC in the distance between the particles!! 
                       if ((dx>0) && (dx>(box/2)))  dx=dx-box;
					        if ((dx<0) && (dx< -(box/2))) dx=dx+box;	
                        if ((dy>0) && (dy>(box/2)))  dy=dy-box;
					        if ((dy<0) && (dy< -(box/2))) dy=dy+box;	
                        
                       
                    //     dx -= box * round(dx / box);
                    //     dy -= box * round(dy / box);

                    // if (i==4) {
                    //     std::cout << "cellX = " << cellX << ", "<< "cell Y = " <<cellY << std::endl;
                    //     std::cout << j << " is neighbor in  cellX = " << cellX << ", "<< "cell Y = " <<cellY<<std::endl;
                    // }

                       

                        double distanceSquared = dx * dx + dy * dy;
                       
                        // if (distanceSquared < CUTOFF * CUTOFF) { 

                        //     // ---- FORCE ---- // 
                            
                        //     double force=0.0;
                            

                        //     double distance = sqrt(distanceSquared); // LJ potential 
                        //     double sigma_12 = pow(sigma,12);

                        //     double sigma12_over_r_11 = sigma_12* pow(distance, 11);
                        //      //std::cout << "sigma/r= " << sigma_over_r_12 << "\n";
    
                        //     double forceMagnitude = 48 * epsilon * (sigma12_over_r_11);     
                        //     if (distance < pow(2,1/6)*sigma) {force = forceMagnitude + epsilon;}
                        //     else {force = 0.0;}
                            

                        //     double fx = force * dx / distance;
                        //     double fy = force * dy / distance;
                            
                            



                            
                        //     pi.fx += fx;
                        //     pi.fy += fy;
                        //     pj.fx -= fx; 
                        //     pj.fy -= fy;
                            
                        //     // std::cout << "Particle " << i << " interacts with Particle " << j << std::endl;
                        // }
                        
                        // double neigh_x = cos(pj.theta)
                        // double neigh_y
                        // std::vector<double> tempArray;
                        if (distanceSquared < R*R){
                            pi.neighbors += 1;
                            double Alignment =  sin(pj.theta - pi.theta);
                            pi.alignment += Alignment; 
                            // pj.alignment -= Alignment; 
                            
                            // neigh_in_region++;
                            order_in_region_x+= cos(pj.theta);
                            order_in_region_y+= sin(pj.theta);
                            
                            pi.neigh_loc.push_back(j);
                            // if ((cos(pj.theta - pi.theta) < -0.9) & (cos(pj.theta - pi.theta) > -1.0)) {
                            //     // std::cout << "particle " << i << " defector found = " << j << std::endl;
                            //     pi.if_defector = 1;
                            //     pi.defector = pj.theta;

                            // } else {pi.if_defector = -1;}

                        }
                        

                        
                            

                        

                    

                    } // loop over neighboring cells
                    order_in_region_x = order_in_region_x/pi.neighbors;
                    order_in_region_y = order_in_region_y/pi.neighbors;
                    pi.avg_ang_region = atan2(order_in_region_y,order_in_region_x);
                    // sqrt(order_in_region_x*order_in_region_x + order_in_region_y+order_in_region_y);
                    pi.op_in_region = sqrt(order_in_region_x*order_in_region_x + order_in_region_y+order_in_region_y);
                    
                    j = linkedList[j];
                } // loop over particles j that intereact with i 

            } // neighboring cells loop y
        } // neighboring cells loop x
	
 } //loop over all particles


}
int time_aux = 0;

// particles_w_high_order.clear();
// particles_w_high_order.shrink_to_fit();
double updatePositions(Particle* particles) {
    //  std::cout << "passed" << std::endl;
    int nbr_poss_defectors=0;
    for (int i = 0; i < N_particles; ++i)
    {
        
        if (particles[i].op_in_region > 0.8) {
            // std::cout << "system ok for chosing defector" << std::endl;
            nbr_poss_defectors ++;
            particles_w_high_order.push_back(i);
        }
    }

int defector;
//  std::cout << "nbr poss defector = "<< nbr_poss_defectors << std::endl;

   if ((nbr_poss_defectors > 0) & (op_now > 0.8)) {
            // std::cout << "defector chosen" << std::endl;
            std::uniform_int_distribution<int> dis_int(0,nbr_poss_defectors-1);
            
            int rand_defector_index = dis_int(gen);
            // std::cout << "here" << std::endl;
            defector = particles_w_high_order[rand_defector_index];

            particles[defector].theta = particles[defector].theta + PI/2;
            // particles[rand_defector_index].theta = particles[rand_defector_index].avg_ang_region + PI/6;
        }

    for (int i = 0; i < N_particles; ++i)
    {
        double update_theta;
        double rnd = dis(gen);
        int nbr_neighbors = particles[i].neighbors;
        std::vector<int> neighbors_of_part_now = particles[i].neigh_loc;
        // for (int j = 0; j < nbr_neighbors-1; j++)
        // {
        //     //  std::cout << "check" << std::endl;
        //     if (neighbors_of_part_now[j] == defector) 
        //     {
        //         break;
        //         double prob_of_defector = dis2(gen);
        //         if (acceptance_rate*Dt >= prob_of_defector) {
        //             update_theta = particles[defector].theta;

        //         } 
        //         else {update_theta = particles[i].alignment*alig_str*Dt/particles[i].neighbors + eta*rnd*sqrt(Dt);}

        //     } else {
        //         update_theta = particles[i].alignment*alig_str*Dt/particles[i].neighbors + eta*rnd*sqrt(Dt);
        //         break;
        //    }
        // }
        std::vector<int>::iterator p;
        p = find(neighbors_of_part_now.begin(),neighbors_of_part_now.end(),defector);
        if (p != neighbors_of_part_now.end())
        {
            // std::cout << "defector found " << std::endl;
           double prob_of_defector = gen();
                if (acceptance_rate*Dt >= prob_of_defector) {
                    // std::cout << "accept rate = " << acceptance_rate*Dt << std::endl;
                    //  std::cout << "prob of defector = " << prob_of_defector << std::endl;
                    // std::cout << "defector activated" << std::endl;
                    update_theta = particles[defector].theta;
                } 
                else {update_theta =  particles[i].theta+ particles[i].alignment*alig_str*Dt/particles[i].neighbors + eta*rnd*sqrt(Dt);}
        } else {
            update_theta =   particles[i].theta+ particles[i].alignment*alig_str*Dt/particles[i].neighbors + eta*rnd*sqrt(Dt);
        }
        
        particles[i].theta = update_theta;
        particles[i].x += particles[i].fx*Dt + Vo*cos(particles[i].theta)*Dt;
        particles[i].y += particles[i].fy*Dt +  Vo*sin(particles[i].theta)*Dt;
        particles[i].aux_posx += particles[i].fx*Dt+ Vo*cos(particles[i].theta)*Dt;
        particles[i].aux_posy += particles[i].fy*Dt+ Vo*sin(particles[i].theta)*Dt;

        particles[i].x = applyPBC(particles[i].x);
        particles[i].y = applyPBC(particles[i].y);
        particles[i].theta = applyPBCangle(particles[i].theta);

        particles[i].cellIndex = getCellIndex(particles[i].x, particles[i].y);
       
    }
    


    




    // double dist=0.0;
    double aux_dist_x=0.0;
    double aux_dist_y=0.0;
    for (int f=0; f <N_particles; f++) {
        // aux_dist_x = particles[f].aux_posx - particles[f].ini_posx;
        // aux_dist_y = particles[f].aux_posy - particles[f].ini_posy;
        // aux_dist_x = aux_dist_x*aux_dist_x;
        // aux_dist_y = aux_dist_y*aux_dist_y;
        // dist += aux_dist_x + aux_dist_y;
        aux_dist_x += cos(particles[f].theta);
        aux_dist_y += sin(particles[f].theta);
       
    }
    aux_dist_x = aux_dist_x/N_particles;
    aux_dist_y = aux_dist_y/N_particles;
    double dist = sqrt(aux_dist_x*aux_dist_x + aux_dist_y*aux_dist_y);
    time_aux++;
    op_now = dist;
    return dist;

}


int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <parameter_file>" << std::endl;
        return 1;
    }

    std::string param_file = argv[1];
    // gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
    // gsl_rng_set(r, time(NULL));

    std::map<std::string, double> params = readParams(param_file);

    // Use the parameters in the simulation
    box = params["box"];
    CUTOFF = params["cutoff"];
    N_particles = static_cast<int>(params["Nparticles"]);
    Vo = params["Vo"];
    eta = params["eta"];
    sigma = params["sigma"];
    R = params["R"];
    epsilon = params["epsilon"];
    alig_str = params["aligstr"];
    T = static_cast<int>(params["T"]);
    Dt = params["Dt"];

    Cells = box / CUTOFF;

    // Generate the dynamic file name
    std::string positions_file_name = createFileName(params, "particle_positions");
    std::string squared_disp_file_name = createFileName(params, "squared_disp");

    Particle* particles = (Particle*)malloc(N_particles * sizeof(Particle)); // use malloc to handle large arrays 
    int* head = (int*)malloc(Cells * Cells * sizeof(int));
    int* linkedList = (int*)malloc(N_particles * sizeof(int));


    if (!particles || !head || !linkedList) {
        std::cerr << "Memory allocation failed" << std::endl;
        return -1;
    }
  
    initializeParticles(particles);
    std::ofstream file(positions_file_name);  // Dynamic output file name for particle positions
    file.close();

    std::ofstream sd;
    sd.open(squared_disp_file_name); 
        // std::ofstream file("particle_positions.dat"); // Create a new file or overwrite if it already exists
        // file.close();
    
    // std::ofstream sd;
    // sd.open("squared_disp.dat");  

    // sd.close();
        
    for (int t = 0; t < T; ++t) {
    //    std::cout << "time= " << t << std::endl << "---------------------" << "\n";
        
  
        buildLinkedList(particles, head, linkedList);


      
   
        simulateInteractions(particles, head, linkedList);



        double store_msd = updatePositions(particles);
        
        sd << (t * Dt) << " " << store_msd << "\n";

            

        // std:: cout << store_msd << "\n";
        // sd << (t*Dt) << " " << store_msd/N_particles << "\n";
        
       if (t%t_XY_Save==0) { 
       savePositions(particles, positions_file_name);}
        
    }

    free(particles);
    free(head);
    free(linkedList);
    sd.close();
    // gsl_rng_free(r);
    return 0;
}
