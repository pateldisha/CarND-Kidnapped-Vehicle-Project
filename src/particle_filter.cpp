/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
		num_particles = 50 ;      //Number of Particles                                                                                              
        default_random_engine gen;                                                                                       
                                                                                                                         
        normal_distribution<double> dist_x(x, std[0]);                                                                   
        normal_distribution<double> dist_y(y, std[1]);                                                                   
        normal_distribution<double> dist_theta(theta, std[2]);                                                           
   
// Initialize particles   
  for(int i = 0; i<num_particles; i++){                                                                                  
    Particle p;                                                                                                          
    p.id = i;                                                                                                            
    p.x = dist_x(gen);                                                                                                   
    p.y = dist_y(gen);                                                                                                   
    p.theta = dist_theta(gen);                                                                                           
    p.weight = 1;                                                                                                      
   weights.push_back(1.0);                                                                                              
    particles.push_back(p);                                                                                              
  }                                                                                                                      
                                                                                                                         
  is_initialized = true;                                                                                                 
  return;      	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	default_random_engine gen;                                                                                       
                                                                                                                         
    normal_distribution<double> dist_x(0, std_pos[0]);                                                               
    normal_distribution<double> dist_y(0, std_pos[1]);                                                               
    normal_distribution<double> dist_theta(0, std_pos[2]);                                                           
                                                                                                                         
  for(int i = 0; i<num_particles; i++){                                                                                                                                                                                                                        

    if(fabs(yaw_rate) > 0.00001){                                                                                        
      particles[i].x = particles[i].x + velocity/yaw_rate * (sin(particles[i].theta + yaw_rate* delta_t)-sin(particles[i].theta));                                               
      particles[i].y = particles[i].y + velocity/yaw_rate * (cos(particles[i].theta)-cos(particles[i].theta + yaw_rate*delta_t));                                               
      particles[i].theta +=yaw_rate*delta_t;                                                                                       
    }else{                                                                                                               
      particles[i].x += velocity * sin(particles[i].theta) * delta_t;                                                            
      particles[i].y += velocity * cos(particles[i].theta) * delta_t;                                                                                   
    }                                                                                                                    
                                                                                                                         
    particles[i].x += dist_x(gen);                                                                                       
    particles[i].y += dist_y(gen);                                                                                       
    particles[i].theta += dist_theta(gen);                                                                               
  }                                                                                                                      

  return;               

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	 weights.clear();                                                                                                 
  for(int i = 0; i<num_particles; i++){  

// Predicted map landmark locations  within the range
    std::vector<LandmarkObs> predicted;                                                                                  

    double p_x = particles[i].x;                                                                                         
    double p_y = particles[i].y;                                                                                         
    double p_theta = particles[i].theta;                                                                                 

    
    particles[i].associations.clear();                                                                                   
    particles[i].sense_x.clear();                                                                                        
    particles[i].sense_y.clear();                                                                      
	// Tranform the observations to map coordinate
    double weight = 1;
	
    for(int j = 0; j<observations.size(); j++){                                                        
	
      double transformed_x = observations[j].x* cos(p_theta) - observations[j].y * sin(p_theta) + p_x; 
      double transformed_y = observations[j].x * sin(p_theta) + observations[j].y * cos(p_theta) + p_y;                                                    
      if(pow(pow(transformed_x-p_x,2)+pow(transformed_y-p_y,2),0.5) > sensor_range) continue;                                        
      particles[i].sense_x.push_back(transformed_x);
      particles[i].sense_y.push_back(transformed_y);                                                                           
      double min_dist = 9999999999;                                                                                     
      int min_id=-1;                                                                                                      
      for(int k = 0; k<map_landmarks.landmark_list.size(); k++){                                                         
                                                                      
        float diff_x = map_landmarks.landmark_list[k].x_f - transformed_x;                                                                                   
        float diff_y = map_landmarks.landmark_list[k].y_f - transformed_y;                                                                                   
        float dist  = pow(pow(diff_x,2)+pow(diff_y,2),0.5);                                                             
        if(dist < min_dist){                                                                                           
          min_dist = dist;                                                                                             
          min_id = k;                                                                                                     
        }                                                                                                                
      }                                                                                                                  
      double ld_x = map_landmarks.landmark_list[min_id].x_f;                                                               
      double ld_y = map_landmarks.landmark_list[min_id].y_f;                                                               

      particles[i].associations.push_back(map_landmarks.landmark_list[min_id].id_i);                                      

      weight = weight * exp(-0.5 * (pow((ld_x - transformed_x) / std_landmark[0],2) + pow((ld_y - transformed_y) / std_landmark[1],2))) / (2*M_PI*std_landmark[0]*std_landmark[1]);                                                                             


    }                                                                                                                    
    particles[i].weight=weight;                                                                                          
    weights.push_back(weight);                                                                                           
  }                                                                                                                      

  return;  
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
default_random_engine gen;                                                                                       
  //discrete_distribution<int> distribution(weights.begin(), weights.end());                                               
  discrete_distribution<int> distribution(weights.begin(), weights.end());                                               
  std::vector<Particle> new_particles;                                                                             

weights.clear();                                                                                                       

  for(int i=0; i < num_particles; i++){                                                                                  
    int j = distribution(gen);                                                                                      
    new_particles.push_back(particles[j]);                                                                    
    weights.push_back(particles[j].weight);                                                                         
  }                                                                                                                      

  particles=new_particles;                                                                                         

  return; 
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
                                                  
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
	    
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
