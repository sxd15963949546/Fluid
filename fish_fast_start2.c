#include "math.h"
#include "unsteady.h"
#include "udf.h"
#include "sg.h"
#include "math.h"
#include "stdio.h"
#include "sg_mem.h"
#include "dynamesh_tools.h"


static real dy = 0;

static real T=0.06;
static real precision = 5e-6; // Precision the calculation
static real A_max=0.265;
static real lamda = 1;
static real num = 1400;
static real t_0=1;
static real t_end=3;
static real pi=M_PI;
static real L_scale = 0.4;
static real xc = 0.50771; // The x-coordinate of the centroid
static real yc = 0; // The y-coordinate of the centroid
#define mass 16

#define orientation -M_PI/6*0
#define orientation_ M_PI/6*0
#define FLUID_DOMAIN 1
#define fish_id 8
double global_zone_total_displacement_x_fish = 0;
double global_zone_velocity_x_fish = 0;
double fish_acceceleration = 0;


double fish_force[3];// the water force on the object
double fish_moment[3];// the moment of the external force on the object with respect to the center position
double fish_center_coordinates[3];// the center coordinates of the fish	
/*
executes at every iteration and is called at the beginning
of every iteration before transport equations are solved.
*/
DEFINE_ADJUST(Adjust_fish1, domain)
{
	if (N_ITER % 15 == 0)//integer number of iterations
	{
		FILE* file_velocity_x;//file_velocity_x
		FILE* file_displacementy_x;

		if ((file_velocity_x = fopen("VelocityX_fish1", "r")) == NULL)//
		{
			Message("Error opening 'VelocityX_fish1 ' for reading...");
		}
		else
		{
			fscanf(file_velocity_x, "%le", &global_zone_velocity_x_fish);//%le 
			fclose(file_velocity_x);
		}

		if ((file_displacementy_x = fopen("TotalDisplacementX_fish1", "r")) == NULL)
		{
			Message("Error opening 'TotalDisplacementX_fish1 ' for reading...");
		}
		else
		{
			fscanf(file_displacementy_x, "%le", &global_zone_total_displacement_x_fish);
			fclose(file_displacementy_x);
		}
	}
}
/*
	The name of the argument that you supplied as the first DEFINE macro argument will become visible and selectable in the cell zone condition dialog boxes
	DEFINE_ZONE_MOTION 
	symbol name 				UDF name.
	real *omega					Pointer to the rotational velocity magnitude, default 0.
	real axis[3]				Rotation axis direction vector, default (0 0 1) and (1 0 0) for 2D axisymmetric cases.
	real origin[3]				Rotation axis origin vector, default (0 0 0).
	real velocity[3]			Translational velocity vector, default (0 0 0).
*/
//In order to use a 3-DOF model in this research, you need to set N3V_S(origin,=,global_zone_total_displacement_x_fish,global_zone_total_displacement_y_fish)
/*
All vectors are specified as 3 dimensional, even for 2D simulations. The third component of the origin and the translational velocity vectors will be ignored in 2D cases. 
For regular 2D cases, rotation is assumed to be around the Z axis. In a 2D axisymmetric case, the rotation is around the X axis.
Hence, for 2D cases any modification to the axis vector inside the UDF will be ignored
*/

DEFINE_ZONE_MOTION(zone_motion, omega, axis, origin, velocity, time, dtime)
{
	Message("Now is zone_motion\n");
	*omega = 0;
	origin[0] = global_zone_total_displacement_x_fish;
	origin[1] = 0;
	velocity[1] = 0;
	velocity[0] = global_zone_velocity_x_fish;
}

DEFINE_GRID_MOTION(point_fish,domain,dt,time,dtime)
{
    Thread *tf = DT_THREAD(dt);
    face_t f;
    Node *v;
	int n, i;	

    double x0, y0;
    double t_past=time-dtime;
    double t_now = time;
    //ArrayList for centerline
    double s_line[1500] = { 0 };
    double x_center_past[1500] = { 0 };
    double y_center_past[1500] = { 0 };
    double slope_past[1500] = { 0 };
    double head_length=0.3;
    double length = 1;
    double seg[1500] =  { 0 };
    int num_head = 0.3*num;
    double x_center_now[1500] = { 0 };
	double y_center_now[1500] = { 0 };
	double slope[1500] = { 0 };
    // ArrayList for outline
	double UP_x_past[1500] = { 0 };
	double UP_y_past[1500] = { 0 };
	double DOWN_x_past[1500] = { 0 };
	double DOWN_y_past[1500] = { 0 };
	double UP_x_now[1500] = { 0 };
	double UP_y_now[1500] = { 0 };
	double DOWN_x_now[1500] = { 0 };
	double DOWN_y_now[1500] = { 0 };
	double Thickness[1500] = { 0 };
    double a_t1;
    double a_x;
    double delta;
    double F,Fx;
    for(i=0;i<num;i++)
       seg[i] = length / num ;
    /*****************************Past  time***********************************/
	/*****************************Centerline***********************************/
	// Use the Newton–Raphson method to update the centerline
    for(i=0; i<=num; i++)
    {
        if(i==0)
        {
            x_center_past[i]=0; 
            y_center_past[i]=0;
        }
        else if(i>0&&i<=0.3*num)
        {
            x_center_past[i] = x_center_past[i-1] + seg[i-1];
            y_center_past[i]=0;
        }
        else
        {
            delta=10000;
            x_center_past[i] = x_center_past[i-1] + seg[i-1];
        
            if(t_past<=t_0)
                a_t1=t_past/t_0-1/(2*pi)*sin(2*pi*t_past/t_0);
            else
            {
                a_t1=(t_end-t_past)/(t_end-t_0)+1/(2*pi)*sin((2*pi*(t_past-t_0)/(t_end-t_0)));
            }
            while (fabs(delta)>precision)
            {
                a_x=A_max*pow(x_center_past[i]-0.3, 2.0);
                y_center_past[i]=a_x*a_t1*sin(2*pi*(x_center_past[i]/lamda-t_past/T));
                
                F = pow(x_center_past[i] - x_center_past[i-1], 2.0) + pow(y_center_past[i] - y_center_past[i-1], 2.0) -pow(seg[i-1], 2.0);
                
                Fx = 2.0 * (x_center_past[i] - x_center_past[i-1]) + 
                2.0 * (y_center_past[i]*((2*A_max*(x_center_past[i]-0.3)*a_t1*sin(2*pi*(x_center_past[i]/lamda-t_past/T)))+a_x*a_t1*(2*pi/lamda)*cos(2*pi*(x_center_past[i]/lamda-t_past/T)))
                - y_center_past[i-1]*((2*A_max*(x_center_past[i-1]-0.3)*a_t1*sin(2*pi*(x_center_past[i-1]/lamda-t_past/T)))+a_x*a_t1*(2*pi/lamda)*cos(2*pi*(x_center_past[i-1]/lamda-t_past/T))));
                
                delta = F / Fx;
                x_center_past[i] = x_center_past[i] - delta;
                /* code */
            } 
            y_center_past[i]=a_x*a_t1*sin(2*pi*(x_center_past[i]/lamda-t_past/T));
        }
    }
	
    /*****************************Outline***********************************/
	for (i = 0; i < num; i++)
	{
        if (i == 0)
		{
			s_line[i] = 0;
		}
        else
        {
            s_line[i] = seg[i] + s_line[i-1];
        }
        slope[i] = (y_center_past[i + 1] - y_center_past[i]) / (x_center_past[i + 1] - x_center_past[i]);
        Thickness[i] = 0.625 * 
                        (0.298222773 * sqrt(s_line[i] ) - 0.127125232 * (s_line[i]) - 0.357907906 * pow(s_line[i] , 2) +
                        0.291984971 * pow(s_line[i] , 3) - 0.105174606 * pow(s_line[i] , 4));

		UP_x_past[i] = x_center_past[i + 1] + Thickness[i] * cos(atan(slope[i]) + M_PI / 2);
		UP_y_past[i] = y_center_past[i + 1] + Thickness[i] * sin(atan(slope[i]) + M_PI / 2);
		DOWN_x_past[i] = x_center_past[i + 1] - Thickness[i] * cos(atan(slope[i]) + M_PI / 2);
		DOWN_y_past[i] = y_center_past[i + 1] - Thickness[i] * sin(atan(slope[i]) + M_PI / 2);
    }

    //UP_x_past[num] = 0;

	//UP_y_past[num] = (a2 * UP_x_past[num] * UP_x_past[num] + a1 * UP_x_past[num] + a0) * sin(k * UP_x_past[num] - t_past * w);

	/****************************Current time**********************************/
	/*****************************Centerline***********************************/
	// Use the Newton–Raphson method to update the centerline
    for(i=0; i<=num; i++)
    {
        if(i==0)
        {
            x_center_now[i]=0; 
            y_center_now[i]=0;
        }
        else if(i>0&&i<=0.3*num)
        {
            x_center_now[i] = x_center_now[i-1] + seg[i-1];
            y_center_now[i]=0;
        }
        else
        {
            delta=10000;
            x_center_now[i] = x_center_now[i-1] + seg[i-1];        
            if(t_now<=t_0)
                a_t1=t_now/t_0-1/(2*pi)*sin(2*pi*t_now/t_0);
            else
            {
                a_t1=(t_end-t_now)/(t_end-t_0)+1/(2*pi)*sin((2*pi*(t_now-t_0)/(t_end-t_0)));
            }
            while (fabs(delta)>precision)
            {
                a_x=A_max*pow(x_center_now[i]-0.3, 2.0);
                y_center_now[i]=a_x*a_t1*sin(2*pi*(x_center_now[i]/lamda-t_now/T));
                
                F = pow(x_center_now[i] - x_center_now[i-1], 2.0) + pow(y_center_now[i] - y_center_now[i-1], 2.0) -pow(seg[i-1], 2.0);
                
                Fx = 2.0 * (x_center_now[i] - x_center_now[i-1]) + 
                2.0 * (y_center_now[i]*((2*A_max*(x_center_now[i]-0.3)*a_t1*sin(2*pi*(x_center_now[i]/lamda-t_now/T)))+a_x*a_t1*(2*pi/lamda)*cos(2*pi*(x_center_now[i]/lamda-t_now/T)))
                - y_center_now[i-1]*((2*A_max*(x_center_now[i-1]-0.3)*a_t1*sin(2*pi*(x_center_now[i-1]/lamda-t_now/T)))+a_x*a_t1*(2*pi/lamda)*cos(2*pi*(x_center_now[i-1]/lamda-t_now/T))));
                
                delta = F / Fx;
                x_center_now[i] = x_center_now[i] - delta;
                y_center_now[i]=a_x*a_t1*sin(2*pi*(x_center_now[i]/lamda-t_now/T));
                /* code */
            } 

        }
    }
    /*****************************Current Time Outline***********************************/
	for (i = 0; i < num; i++)
	{
        if (i == 0)
		{
			s_line[i] = 0;
		}
        else
        {
            s_line[i] = seg[i] + s_line[i-1];
        }
        	slope[i] = (y_center_now[i + 1] - y_center_now[i]) / (x_center_now[i + 1] - x_center_now[i]);
		    Thickness[i] = 0.625 * 
			(0.298222773 * sqrt(s_line[i] ) - 0.127125232 * (s_line[i] ) - 0.357907906 * pow(s_line[i] , 2) +
				0.291984971 * pow(s_line[i] , 3) - 0.105174606 * pow(s_line[i] , 4));

		UP_x_now[i] = x_center_now[i + 1] + Thickness[i] * cos(atan(slope[i]) + M_PI / 2);
		UP_y_now[i] = y_center_now[i + 1] + Thickness[i] * sin(atan(slope[i]) + M_PI / 2);
		DOWN_x_now[i] = x_center_now[i + 1] - Thickness[i] * cos(atan(slope[i]) + M_PI / 2);
		DOWN_y_now[i] = y_center_now[i + 1] - Thickness[i] * sin(atan(slope[i]) + M_PI / 2);
    }
    // FILE* file_x_up_data, file_y_up_data;
    // FILE* file_x_down_data, file_y_down_data;
    // file_x_up_data = fopen("X_UP_DATA.txt","w");
    // file_y_up_data = fopen("Y_UP_DATA.txt","w");
    // file_x_down_data = fopen("X_DOWN_DATA.txt","w");
    // file_y_down_data = fopen("Y_DOWN_DATA.txt","w");
    //UP_x_now[num] = 0;
	//UP_y_now[num] = (a2 * UP_x_now[num] * UP_x_now[num] + a1 * UP_x_now[num] + a0) * sin(k * UP_x_now[num] - t_now * w);
	/*****************************Update***********************************/
    /* set deforming flag on adjacent cell zone */
    SET_DEFORMING_THREAD_FLAG(THREAD_T0(tf));
    Message("Now is mesh_motion\n");
    begin_f_loop(f,tf)
    {
        f_node_loop(f,tf,n)
        {

            v = F_NODE(f,tf,n);
            /* update node */
            if( NODE_POS_NEED_UPDATE(v))
            {
                NODE_POS_UPDATED(v);
				double xt, yt, x0, y0;
				double Distance_up = 1e+03;
				double Distance_down = 1e+03;
				double Distance;
				double x_diff, y_diff;
				int sign_up, sign_down;

				//Change the coordinates to horizontal coordinates
				xt = NODE_X(v) / L_scale - global_zone_total_displacement_x_fish/L_scale ; 
                yt = NODE_Y(v) / L_scale - dy/L_scale ;
				x0 = xt; y0 = yt;
				for (i = 0; i < num; i++)
				{
					x_diff = xt - UP_x_past[i];
					y_diff = yt - UP_y_past[i];
					Distance = sqrt((x_diff * x_diff) + (y_diff * y_diff));
					if (Distance <= Distance_up)
					{
						sign_up = i;
						Distance_up = Distance;
					}
				}
				
				for (i = 0; i < num; i++)
				{
					x_diff = xt - DOWN_x_past[i];
					y_diff = yt - DOWN_y_past[i];
					Distance = sqrt((x_diff * x_diff) + (y_diff * y_diff));
					if (Distance <= Distance_down)
					{
						sign_down = i;
						Distance_down = Distance;
					}
				}

				if (Distance_up < Distance_down)
				{
					NODE_X(v) = NODE_X(v) + (UP_x_now[sign_up] - UP_x_past[sign_up]) * L_scale;
					NODE_Y(v) = NODE_Y(v) + (UP_y_now[sign_up] - UP_y_past[sign_up]) * L_scale;
                   
				}
				else
				{
					NODE_X(v) = NODE_X(v) + (DOWN_x_now[sign_down] - DOWN_x_past[sign_down]) * L_scale;
					NODE_Y(v) = NODE_Y(v) + (DOWN_y_now[sign_down] - DOWN_y_past[sign_down]) * L_scale;
				}
            }            
        }

    }
	end_f_loop(t,tf);    
}
/*In this part, I call the Compute_force_and_moment macro which referenced the follow website 
	https://zhuanlan.zhihu.com/p/33577479
		When I used this macro in parallel mode, I found that there is no need to calculate Force and Moment 
	in different nodes, which means you can get the same results in different nodes.
		The logical variables TRUE and FALSE at the end of the macro do not indicate whether to call the macro. 
	Instead, it indicates whether to run the macro in the main node host to get the force and moment calculation 
	results also in the main node.
*/
DEFINE_EXECUTE_AT_END(calculate_velocity_fish)
{
/*Different variebles are needed on different nodes*/
/*Variebles that defined in host*/
    double dtime = CURRENT_TIMESTEP;
	double time = CURRENT_TIME;

	//The variables for host and node communication must be global variables

/*Variebles that defined in nodes*/
#if RP_NODE
	Domain* d = Get_Domain(FLUID_DOMAIN);
    Thread* f_thread;
#endif

/*Define file type variable*/
#if RP_HOST
	FILE* file_displacementy_x;
    FILE* file_velocity_x;
	FILE* file_acceceleration_x;
    FILE* file_results_fish;
#endif

#if RP_NODE/*Nodes are used for calculations such as element data*/
	if (!Data_Valid_P())
		return;
    f_thread = Lookup_Thread(d,fish_id);
    Compute_Force_And_Moment(d, f_thread, fish_center_coordinates, fish_force, fish_moment,0);
 //   fish_force[0] = fish_force[0];

	real global_zone_velocity_x_fish_last = global_zone_velocity_x_fish;
	real fish_acceceleration_last = fish_acceceleration;
    fish_acceceleration = (fish_force[0]) / mass;
    global_zone_velocity_x_fish = global_zone_velocity_x_fish + (fish_acceceleration_last+fish_acceceleration) * dtime/2;
    global_zone_total_displacement_x_fish = global_zone_total_displacement_x_fish + (global_zone_velocity_x_fish + global_zone_velocity_x_fish_last) * dtime/2;

#endif
    node_to_host_real_7(fish_force[0],fish_force[1],
						fish_moment[0],fish_moment[1],fish_moment[2],
						global_zone_velocity_x_fish,global_zone_total_displacement_x_fish);

/*Host is used to write or read documents with outside*/
#if RP_HOST
	file_results_fish = fopen("udf-results-fish.dat", "a+");
	fprintf(file_results_fish, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n", time, global_zone_velocity_x_fish,global_zone_total_displacement_x_fish,
			fish_force[0], fish_force[1], fish_moment[0], fish_moment[1], fish_moment[2]);
	fclose(file_results_fish);

	file_velocity_x = fopen("VelocityX_fish1", "w");
	fprintf(file_velocity_x, "%le", global_zone_velocity_x_fish);
	fclose(file_velocity_x);

	file_displacementy_x = fopen("TotalDisplacementX_fish1", "w");
	fprintf(file_displacementy_x, "%le", global_zone_total_displacement_x_fish);
	fclose(file_displacementy_x);
#endif
}


/*
DEFINE_ON_DEMAND(assignment)
{
	real x_date[4] = { 0,0,0,0 };

	dx1 = x_date[3];
	x = x_date[0];
	ve = x_date[1];
	a = x_date[2];

	Message("dx=%f,x=%f,vx=%f,ax=%f \n",dx1,x,ve,a);
} 
*/