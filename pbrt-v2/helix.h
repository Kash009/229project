class Helix{
  public:
  	Helix(float thickness, unsigned int ridge, float n_twirls, float twist, float r_start, float r_end);
  	void Mesh(unsigned int nsamples);	
  private:
	float d(float t);
	float r(float t);
	float t_end;//is n_twirls + 1.5(half a twirl for beginning, half for end)
	float thickness;
	float r_start;
	float r_increment;
	//for cream/softserve/etc
	float type;
	//for sampling
	int sample_rate_cs;
	int sample_rate_cv;
	//RandomSampler
}

//distance to z-axis
inline float Helix::d(float t){
  if(t <1)return r_start;	
  if(t <t_end-0.5f) return r_start +s*r_increment;
  return (r_start+r_increment)*powf(2*(t_end-t),2);
}

//height (beginning at 0)
float Helix::h(float t){
  //first layer
  if(t <1)return (1+t)*thickness;
  //middle layers
  if(t<t_end-0.5f) return h(t-1)+sqrt(thickness*thickness - powf((d(t)-d(t-1)),2));
  return -(t_end-t);
}

inline float Helix::twist(float t){

}
