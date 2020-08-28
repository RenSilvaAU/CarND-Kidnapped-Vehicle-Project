/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"


using std::string;
using std::vector;
using std::normal_distribution;
using std::default_random_engine;

static default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO - DONE: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO - DONE: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 130;  // TODO - DONE: Set the number of particles

	Particle p;
  p.weight = 1; 

	// start up with "noisy" normal distribution
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

	for(int i=0; i < num_particles; i++) {

    // define particle
    p.id=i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
 
    // add to vector
		particles.push_back(p);
		weights.push_back(1);

	}

	is_initialized = true;
}



void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO - DONE: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  // for every particle
	for(int i=0; i < num_particles; i++) {

		if (abs(yaw_rate) > NEAR_ZERO) {

      // predict
			particles[i].x = particles[i].x + velocity * (sin(particles[i].theta + yaw_rate * delta_t) 
                      - sin(particles[i].theta)) / yaw_rate;

			particles[i].y = particles[i].y + velocity * (cos(particles[i].theta) 
                      - cos(particles[i].theta + yaw_rate * delta_t)) / yaw_rate;

			particles[i].theta = particles[i].theta + yaw_rate * delta_t;

		} else {

      // if yaw_rate tends to zero, simplify the equations
			particles[i].x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			particles[i].y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
      // particles[i].theta is unchanged

		}

    // now add some noise
		normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
		normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
		normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

	}
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO - DONE: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */


	for (int i = 0; i < observations.size(); i++) {

		double mininum_distance = 0;

		for (int j = 0; j < predicted.size(); j++) {

			double distance = dist(observations[i].x, observations[i].y, 
                             predicted[j].x,    predicted[j].y);

			if (j == 0) 
				mininum_distance = distance;

			else 

        // find the landmark obs with minimum distance
				if (distance < mininum_distance) {

					mininum_distance = distance;

          if (predicted[j].id > 0)
					  observations[i].id = predicted[j].id-1;

				}
			
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   *  TODO - DONE: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

	double gaussian_normalised = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);

	for (int i = 0; i < num_particles; i++) {

    // predicted landmarks
		vector<LandmarkObs> predicted_landmarks;
		
    for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {

      LandmarkObs landmark_obs;

      landmark_obs.id = map_landmarks.landmark_list[j].id_i;
      landmark_obs.x = map_landmarks.landmark_list[j].x_f;
      landmark_obs.y = map_landmarks.landmark_list[j].y_f;

      predicted_landmarks.push_back(landmark_obs);

    }
		
    // observations map
		vector<LandmarkObs> obs_map;
		
    for (int k = 0; k < observations.size(); k++) {

      LandmarkObs trans_lmo;
      
      trans_lmo.id = observations[k].id;
      trans_lmo.x = particles[i].x + observations[k].x * cos(particles[i].theta) 
                  - observations[k].y * sin(particles[i].theta);
      trans_lmo.y = particles[i].y + observations[k].x * sin(particles[i].theta) + observations[k].y * cos(particles[i].theta);
      obs_map.push_back(trans_lmo);
      
    }
		
    // now associate them
		dataAssociation(predicted_landmarks, obs_map);

		
    double prob = 1.0;

    // let us update the weights
    for( int k = 0; k < obs_map.size(); k++) {

      if (obs_map[k].id > predicted_landmarks.size() ) 
        continue;

      LandmarkObs predicted = predicted_landmarks[obs_map[k].id];

      double distance = dist(predicted.x, predicted.y, particles[i].x, particles[i].y);

      if (distance > sensor_range)
        continue;

      double mu_x = predicted.x, mu_y = predicted.y;

      double exponent = pow( obs_map[k].x - mu_x, 2 ) / (2 * pow(std_landmark[0], 2) ) 
                      + pow( obs_map[k].y - mu_y, 2 ) / (2 * pow(std_landmark[1], 2) );

      prob *= gaussian_normalised * exp( -1 * exponent );

    }

    particles[i].weight = prob;
    weights[i] = particles[i].weight;

	}

}

void ParticleFilter::resample() {
  /**
   * TODO - CHECK: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

	std::discrete_distribution<int> discrete_distrib(weights.begin(), weights.end());
	vector<Particle> particles_resample;

	for(int i=0; i < num_particles; i++) {
		particles_resample.push_back( particles.at( discrete_distrib(gen) ) );
	}

	particles = particles_resample;


}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}