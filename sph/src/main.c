#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "particles.h"
#include "acceleration.h"
#include "integrator.h"
#include "initial_conditions.h"
#include "psrand.h"
#include "log.h"
#include "raylib.h"
#include "rcamera.h"
#include "raymath.h"

float plot_rho(particles parts, float h, float R, float rmax, float rho_max, int width, int height, float* rhos, int paused){
    int axis_size = 50;

    int tick_size = 5;

    int samples = 1000;

    float rmin = 0;
    float dx = (rmax - rmin) / ((float)samples - 1.0f);

    int last_x = -100;
    int last_y = -100;

    int nxticks = 7;
    int nyticks = 11;

    float out = 0;
    int zero_found = 0;

    //float rhos[parts.n_particles];
    if (!paused){
        #pragma omp parallel for
        for (int i = 0; i < samples; i++){
            float r = rmin + ((float)i) * dx;

            float total_rho = 0;

            vec3 vec_r = v3(r,0,0);
            float rho = getRho(vec_r,parts,h);

            rhos[i*3] = rho;

            if (rho > rho_max){
                rho_max = rho;
            }

            total_rho += rho;

            //if (((rho/rhos[0]) <= 0.2) && (!zero_found)){
            //    zero_found = 1;
            //    out = r;
            //}

            vec_r = v3(0,r,0);
            rho = getRho(vec_r,parts,h);

            rhos[i*3 + 1] = rho;

            if (rho > rho_max){
                rho_max = rho;
            }

            total_rho += rho;

            //if (((rho/rhos[0]) <= 0.2) && (!zero_found)){
            //    zero_found = 1;
            //    out = r;
            //}

            vec_r = v3(0,0,r);
            rho = getRho(vec_r,parts,h);

            rhos[i*3 + 2] = rho;

            if (rho > rho_max){
                rho_max = rho;
            }

            total_rho += rho;

            if ((fabs(total_rho - 0.0001f) < 0.00005f) && (!zero_found)){
                zero_found = 1;
                out = r;
            }

            //if (zero_found && (total_rho >= 0.0000001f)){
            //    zero_found = 0;
            //    out = r;
            //}


        }
    }

    Vector2 points[3][samples];
    for (int i = 0; i < samples; i++){
        float r = rmin + ((float)i) * dx;
        int render_r = ((r - rmin) / (rmax - rmin)) * (width - axis_size*2) + axis_size;

        points[0][i].x = render_r;
        points[1][i].x = render_r;
        points[2][i].x = render_r;

        //vec3 vec_r = v3(r,0,0);
        float rho = rhos[i*3];
        int render_rho = height - ((rho / rho_max) * (height - axis_size*2) + axis_size);
        points[0][i].y = render_rho;

        rho = rhos[i*3 + 1];
        render_rho = height - ((rho / rho_max) * (height - axis_size*2) + axis_size);
        points[1][i].y = render_rho;

        rho = rhos[i*3 + 2];
        render_rho = height - ((rho / rho_max) * (height - axis_size*2) + axis_size);
        points[2][i].y = render_rho;
    }

    DrawLineStrip(points[2],samples,YELLOW);
    DrawLineStrip(points[0],samples,GREEN);
    DrawLineStrip(points[1],samples,RED);

    DrawLine(axis_size,axis_size,axis_size,height - axis_size,WHITE);
    DrawLine(axis_size,height - axis_size,width - axis_size,height - axis_size,WHITE);

    dx = (rmax - rmin) / ((float)nxticks - 1.0f);
    for (int i = 0; i < nxticks; i++){
        float r = dx * ((float)i) + rmin;
        int render_r = ((r - rmin) / (rmax - rmin)) * (width - axis_size*2) + axis_size;
        //DrawCircle(render_r,axis_size,5,RED);
        DrawLine(render_r,height - axis_size,render_r,(height - axis_size) + tick_size,WHITE);
        char x_text[50];
        sprintf(x_text,"%.2f",r/R);
        DrawText(x_text,render_r,(height - axis_size) + tick_size,10,WHITE);
    }

    float dy = rho_max / ((float)nyticks - 1.0f);
    for (int i = 0; i < nyticks; i++){
        float rho = dy * ((float)i);
        int render_rho = height - ((rho / rho_max) * (height - axis_size*2) + axis_size);
        DrawLine(axis_size - tick_size,render_rho,axis_size,render_rho,WHITE);
        char y_text[50];
        sprintf(y_text,"%.2f",rho);
        int string_width = ((float)MeasureText(y_text,10)) * 1.2f;
        int string_height = MeasureTextEx(GetFontDefault(),y_text,10,1).y;
        DrawText(y_text,(axis_size - tick_size) - string_width,render_rho - (string_height/2),10,WHITE);
    }

    return out;

}


int shader_loc;
Shader shader;

void draw_particles(particles parts, Mesh sphere_mesh, Material mat){
    Vector3* pos = (Vector3*)parts.pos;
    vec3* vel = parts.vel;
    int n_particles = parts.n_particles;

    float R = 0;
    float speeds[n_particles];
    float rs[n_particles];
    float max_speed = 0.01;
    for (int i = 0; i < n_particles; i++){
        speeds[i] = v3len(vel[i]);
        float r = v3len(parts.pos[i]);
        rs[i] = r;
        if (r > R){
            R = r;
        }
        if (speeds[i] > max_speed){
            max_speed = speeds[i];
        }
    }

    SetShaderValue(shader,shader_loc,&R,SHADER_UNIFORM_FLOAT);

    Matrix transforms[n_particles];

    for (int i = 0; i < n_particles; i++){
        float color = rs[i]/R;//speeds[i] / max_speed;
        //if (speed > 1.0f){
        //    LOG("SPEED TOO HIGH!!");
        //}

        Color col = WHITE;
        //col.r = ((float)col.r) * speed;
        col.g = ((float)col.g) * (1.0f-color);
        col.b = 0;//((float)col.b) * (1.0f-color);
        col.a = 200;

        Matrix translation = MatrixTranslate(pos[i].x,pos[i].y,pos[i].z);
        Vector3 axis = Vector3One();
        float angle = 0;
        Matrix rotation = MatrixRotate(axis, angle);

        transforms[i] = MatrixMultiply(rotation, translation);
        DrawMesh(sphere_mesh,mat,transforms[i]);

        //transforms[i] =

        //DrawSphere(pos[i],r,col);
    }

    Color col = WHITE;
    col.a = 20;
    DrawSphereWires(Vector3Zero(),R,100,100,col);

}

float camera_pitch = 0;

void RotateCamera(Camera *camera, float speed)
{
    // Orbital can just orbit
    Matrix rotation = MatrixRotate(camera->up, (IsKeyDown(KEY_RIGHT) - IsKeyDown(KEY_LEFT))*speed*GetFrameTime());
    Vector3 view = Vector3Subtract(camera->position, camera->target);
    view = Vector3Transform(view, rotation);
    camera->position = Vector3Add(camera->target, view);

    CameraPitch(camera, (IsKeyDown(KEY_DOWN) - IsKeyDown(KEY_UP))*speed*GetFrameTime(), 1, 1, 0);

    Vector3 forward = Vector3Scale(Vector3Normalize(Vector3Subtract(camera->target,camera->position)),GetMouseWheelMove());
    camera->position = Vector3Add(camera->position,forward);
}

int main(){
    const int screenWidth = 1000;
    const int screenHeight = 1000;
    InitWindow(screenWidth, screenHeight, "sph");

    const int plotWidth = 400;
    const int plotHeight = 200;

    RenderTexture2D plot_tex = LoadRenderTexture(plotWidth,plotHeight);

    Camera camera = { 0 };
    camera.position = (Vector3){ 2.0f, 2.0f, 2.0f };    // Camera position
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera projection type

    set_seed(21082001);

    float M = 20;
    int N = 20000;
    float R = 1;

    float h = 0.01f;//0.04 / sqrtf((float)N / 1000.0f);
    float k = 0.1f;
    float n = 5;
    //float lambda = 2.01;
    float nu = 1;
    float dt = 0.01;
    float mass = M/(float)N;
    float t = 0;

    float lambda = 2.0f*k*(1.0f+n)*powf(M_PI,(-3.0f/(2.0f*n))) * powf((M*tgammaf(5.0f/2.0f+n)/(R*R*R)/tgammaf(1.0f+n)),(1.0f/n)) / (R*R);
    //printf("lambda = %g\n",lambda);

    float timestep_controller = -0.1;

    particles parts = particles_alloc(N,mass);
    initialize_random_sphere(parts,R * 3);
    //initialize_random(parts,-0.5,0.5,-0.5,0.5,-0.5,0.5);

    leapfrog_init(parts,h,k,n,lambda,nu);

    Mesh sphere_mesh = GenMeshSphere(((2.0f/(float)N)/0.0025f) * 0.01f,10,10);

    shader = LoadShader("/Users/humzaqureshi/GitHub/stars-project/sph/src/shaders/particles.vs","/Users/humzaqureshi/GitHub/stars-project/sph/src/shaders/particles.fs");

    Material mat = LoadMaterialDefault();
    mat.maps[MATERIAL_MAP_DIFFUSE].color = BLUE;
    mat.shader = shader;
    shader_loc = GetShaderLocation(shader,"R");

    //particles_debug(parts);

    int paused = 1;
    float rhos[1000 * 3];
    int computed = 0;

    float sphr = R;

    while (!WindowShouldClose()){

        RotateCamera(&camera,0.5);

        if (IsKeyPressed(KEY_SPACE)){
            paused = !paused;
        }

        BeginDrawing();

        BeginTextureMode(plot_tex);
        Color clear_back = WHITE;
        clear_back.a = 0;
        ClearBackground(clear_back);
        float tmp = plot_rho(parts,h,R,R * 1.5f,5,plotWidth,plotHeight,rhos,paused && computed);
        if (!(paused && computed)){
            sphr = tmp;
        }
        computed = 1;
        EndTextureMode();

        ClearBackground(BLACK);

        BeginMode3D(camera);

            draw_particles(parts,sphere_mesh,mat);

            Color star_col = YELLOW;
            star_col.a = 10;

            DrawSphereWires(Vector3Zero(),sphr,100,100,star_col);

            DrawGrid(10,1);

        EndMode3D();

        //DrawTexture(plot_tex.texture,GetRenderWidth() - plotWidth,0,WHITE);
        DrawTextureRec(plot_tex.texture, (Rectangle){ 0, 0, (float)plot_tex.texture.width, (float)-plot_tex.texture.height }, (Vector2){ GetRenderWidth() - plotWidth, 0 }, WHITE);

        DrawFPS(10,10);

        char time_text[50];
        sprintf(time_text,"time = %.2f",t);
        DrawText(time_text,10,50,10,RED);

        if (paused){
            DrawText("Paused",10,100,10,RED);
        } else {
            DrawText("Playing",10,100,10,GREEN);
        }

        EndDrawing();

        if (!paused){
            timestep_controller += GetFrameTime();
            if (timestep_controller > dt){
                leapfrog_timestep(parts,dt,h,k,n,lambda,nu);
                t += dt;
                timestep_controller = 0;
            }
        }
    }

    particles_destroy(parts);
    UnloadMesh(sphere_mesh);
    UnloadShader(shader);

    CloseWindow();
    return 0;
}