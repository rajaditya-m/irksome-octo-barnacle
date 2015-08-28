#include "RenderObject.h"

RenderObject::RenderObject(SceneObjectId id,bool rotInv)
{
	sceneId_ = id;
	rotationInvariant_ = rotInv;
}

RenderObject::~RenderObject()
{
}

void RenderObject::render(int frameId)
{
	if(render_)
	{
	 
		switch(renderMode_)
		{
		case WIREFRAME : mesh_->draw_wireframe(frameId); break;
		case SHADING : mesh_->draw_mesh(frameId);  break; 
		default : std::cout << "[ERROR] The selected mode of rendering is not supported for this object.\n"; break;
		}
	}
}