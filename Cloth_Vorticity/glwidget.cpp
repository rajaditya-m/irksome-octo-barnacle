#include "glwidget.h"
#include "OpenGLHelperFuncs.cpp"

//Using stuff in Clothing_3 before

GLWidget::GLWidget(QWidget *parent)
	: QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
	//Initialize scene
	scene_ = new Scene();

	//body and cloth data - renderable objects 
	/*body_information_ = new Body_Data("Sphere"
		,"C:/Users/Rajaditya/Copy/Models/Drop_Cloth_Test/Collision_Obj_vertex.bin",
		"C:/Users/Rajaditya/Copy/Models/Drop_Cloth_Test/Collision_Obj_mesh.bin",
		"C:/Users/Rajaditya/Copy/Models/Drop_Cloth_Test/Collision_Obj_textures.bin",
		"body_material.xml");*/
	body_information_ = new Body_Data("Sphere"
		,"D://Copy/Cloth_OBJS/support.obj",
		"body_material.xml");
	body_information_->setRenderable(true);
	body_information_->setRenderMode(SHADING);
	scene_->addRenderObject(body_information_);
	std::cout << "[INFO] Body Information read successfully.\n";
	/*cloth_information_ = new Cloth_Data("Drop_Cloth","C:/Users/Rajaditya/Copy/Models/Drop_Cloth_Test/Cloth_Sheet_vertex.bin",
		"C:/Users/Rajaditya/Copy/Models/Drop_Cloth_Test/Cloth_Sheet_mesh.bin",
		"C:/Users/Rajaditya/Copy/Models/Drop_Cloth_Test/Cloth_Sheet_textures.bin",
		"cloth_physical_param.xml",
		"cloth_material.xml");*/
	cloth_information_ = new Cloth_Data("Drop_Cloth",
		"D://Copy/Cloth_OBJS/cloth2.obj",
		"cloth_physical_param.xml",
		"cloth_material.xml");
	cloth_information_->setRenderable(true);
	cloth_information_->setRenderMode(SHADING);
	scene_->addRenderObject(cloth_information_);
	std::cout << "[INFO] Cloth Information read successfully.\n";

	//Solver simu lation engine and collision engine data
	fem_solver_ = new FEM_Solver();
	collisionEngine_ = new CollisionEngine();
	sim_engine_ = new SimulationEngine(cloth_information_,fem_solver_,body_information_);



	//Animaton timer
	animation_timer = new QTimer(this);
	connect(animation_timer, SIGNAL(timeout()), this, SLOT(nextFrame()));

	//Generate the voxel data 
	float cloth_avg_len = cloth_information_->getMesh()->get_average_edge_length();
	float body_avg_len = cloth_information_->getMesh()->get_average_edge_length();
	float grid_resolution = std::min(cloth_avg_len,body_avg_len)*5.0f;
	//std::cout << "[INFO] Grid resolution is " << grid_resolution << "\n";

	Pair_Vec3d body_super_bbox = body_information_->getMesh()->get_super_bounding_box();
	Eigen::Vector3d min_s_box = body_super_bbox.first;
	float thickness = cloth_information_->get_property_obj()->get_thickness();
	Eigen::Vector3d adjustment(thickness*10.0f,thickness*10.0f,thickness*10.0f);
	min_s_box -= adjustment;

	cloth_information_->getMesh()->get_grid_obj()->set_grid_parameters(min_s_box,grid_resolution);
	body_information_->getMesh()->get_grid_obj()->set_grid_parameters(min_s_box,grid_resolution);

	cloth_information_->getMesh()->update_voxelization();
	//@TODO : Body is not voxelized still now. Possibly add it. 
	//body_information_->getMesh()->voxelize_frame(0);
	//body_information_->getMesh()->voxelize_frame(1);

}


GLWidget::~GLWidget()
{
}

QSize GLWidget::minimumSizeHint() const
{
		return QSize(50, 50);
}

QSize GLWidget::sizeHint() const
{
		return QSize(400, 400);
}

void GLWidget::initializeGL()
{
	scene_->setupScene();
}

void GLWidget::paintGL()
{
	scene_->renderScene();
}

void GLWidget::resizeGL(int width,int height)
{
	scene_->resizeScene(width,height);
}

static void qNormalizeAngle(int &angle)
{
	while(angle < 0)
	{
		angle += 360*16;
	}
	while(angle > 360*16)
	{
		angle -= 360*16;
	}
}

void GLWidget::setRenderingFrameNumber(double frame_number)
{
	int frame_nos = static_cast<int>(frame_number);
	scene_->setCurrentRenderFrame(frame_nos);
	int renderFrame = frame_nos;
	int cloth_frame_available = (cloth_information_->getMesh()->get_number_frames())-1;
	if(renderFrame>cloth_frame_available)
	{
		int diff = renderFrame-cloth_frame_available;
		for(int num = 0; num < diff ; num++)
			sim_engine_->generate_next_frame();
	}
	emit frameNumberUpdated(frame_nos);
	updateGL();
}

void GLWidget::setRenderingFrameNumber(int frame_number)
{
	int renderFrame = frame_number;
	scene_->setCurrentRenderFrame(renderFrame);
	int cloth_frame_available = (cloth_information_->getMesh()->get_number_frames())-1;
	if(renderFrame>cloth_frame_available)
	{
		int diff = renderFrame-cloth_frame_available;
		for(int num = 0; num < diff ; num++)
			sim_engine_->generate_next_frame();
	}
	double frame_nos = static_cast<double>(frame_number);
	emit frameNumberUpdated(frame_nos);
	updateGL();
}

void GLWidget::startAnimation()
{
	animation_timer->start(30);
}

void GLWidget::pauseAnimation()
{
	animation_timer->stop();
}

void GLWidget::resetAnimation()
{
	animation_timer->stop();
	scene_->setCurrentRenderFrame(0);
	emit animFrameNumberUpdated(0);
	emit frameNumberUpdated(0.0); 
	cloth_information_->resetParameters();
	fem_solver_->resetParameters();
	sim_engine_->populatePerVertexBuffer();
	updateGL();
}

//Signals for the 'Render' Menu options 

void GLWidget::bodyDisplayToggled(bool checked)
{
	body_information_->setRenderable(checked);
	updateGL();
}

void GLWidget::clothDisplayToggled(bool checked)
{
	cloth_information_->setRenderable(checked);
	updateGL();
}

void GLWidget::renderGroundToggled(bool checked) 
{
	scene_->setRenderGround(checked);
	updateGL();
}

void GLWidget::renderAxesToggled(bool checked)
{
	scene_->setRenderAxes(checked);
	updateGL();
}

void GLWidget::setBodyRenderInShading(bool checked)
{
	body_information_->setRenderMode(SHADING);
	updateGL();
}

void GLWidget::setBodyRenderInWireframe(bool checked)
{
	body_information_->setRenderMode(WIREFRAME);
	updateGL();
}

void GLWidget::setClothRenderInShading(bool checked)
{
	cloth_information_->setRenderMode(SHADING);
	updateGL();
}

void GLWidget::setClothRenderInWireframe(bool checked)
{
	cloth_information_->setRenderMode(WIREFRAME);
	updateGL();
}

void GLWidget::setClothRenderInHeatmapVelocity(bool checked)
{
	cloth_information_->setRenderMode(HMAP_VELOCITY);
	updateGL();
}

void GLWidget::setClothRenderInHeatmapAcceleration(bool checked)
{
	cloth_information_->setRenderMode(HMAP_ACCLRN);
	sim_engine_->populatePerVertexBuffer();
	updateGL();
}

void GLWidget::setClothRenderInHeatmapBendingForce(bool checked)
{
	cloth_information_->setRenderMode(HMAP_BENDING_FORCE);
	sim_engine_->populatePerVertexBuffer();
	updateGL();
}

void GLWidget::setClothRenderInHeatmapDampingForce(bool checked)
{
	cloth_information_->setRenderMode(HMAP_DAMPING_FORCE);
	sim_engine_->populatePerVertexBuffer();
	updateGL();
}

void GLWidget::setClothRenderInHeatmapShearForce(bool checked)
{
	cloth_information_->setRenderMode(HMAP_SHEAR_FORCE);
	sim_engine_->populatePerVertexBuffer();
	updateGL();
}

void GLWidget::nextFrame()
{
	int renderFrame = scene_->getCurrentRenderFrame();
	renderFrame++;
	if(renderFrame==301)
		renderFrame = 0;
	scene_->setCurrentRenderFrame(renderFrame);
	int cloth_frame_available = (cloth_information_->getMesh()->get_number_frames())-1;
	if(renderFrame>cloth_frame_available)
	{
		int diff = renderFrame-cloth_frame_available;
		for(int num = 0; num < diff ; num++)
			sim_engine_->generate_next_frame();
	}
	emit animFrameNumberUpdated(renderFrame);
	emit frameNumberUpdated((double)renderFrame);
	updateGL();
}

void GLWidget::mousePressEvent(QMouseEvent* event)
{
	scene_->setLastPosition(event->pos());
}

void GLWidget::mouseMoveEvent(QMouseEvent* event)
{
	QPoint cur_pos = event->pos();
	int xRot = scene_->getXRot();
	int yRot = scene_->getYRot();
	float scaling = scene_->getScaling();

	QPoint last_pos = scene_->getLastPosition();
	QPoint diff = cur_pos - last_pos;
	
	if(event->buttons() & Qt::RightButton)
	{
		int x_angle = xRot + 8*diff.y();
		qNormalizeAngle(x_angle);
		if(x_angle != xRot)
		{
			scene_->setXRot(x_angle);
		}
		int y_angle = yRot + 8*diff.x();
		qNormalizeAngle(y_angle);
		if(y_angle != yRot)
		{
			scene_->setYRot(y_angle);
		}
	}
	else if(event->buttons() & Qt::MiddleButton)
	{
		float old_size = scaling;
		scaling *= (1 + (diff.y())/60.0);
		scene_->setScaling(scaling);
		if (scaling <0)
		{
			scene_->setScaling(old_size);
		}
	}
	else if(event->buttons() & Qt::LeftButton) {
		float depth = getPixelDepth(cur_pos.x(),cur_pos.y());
	}
	updateGL();
	scene_->setLastPosition(cur_pos);
}

