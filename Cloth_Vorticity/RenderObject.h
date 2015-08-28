#pragma once

#include "triMesh.h"
#include "mesh_material.h"
#include "global_typedefs.h"

//A renderable object has the following things in it :
// A mesh and a material xml file  


class RenderObject 
{
public:
	RenderObject(SceneObjectId id,bool rotInv);
	~RenderObject();

	virtual void render(int frameId);

	void setMesh(TriMesh* mesh)					{ mesh_ = mesh;					}
	TriMesh* getMesh() const					{ return mesh_;					}

	void setRenderable(bool render)				{ render_ = render;				}
	bool isRenderable() const					{ return render_;				}

	void setRenderMode(RenderMode rendermode)	{ renderMode_ = rendermode;		}
	RenderMode getRenderMode() const			{ return renderMode_;			}

	SceneObjectId getSceneObjectId()			{ return sceneId_;				}

	bool isRotationInvariant() const			{ return rotationInvariant_;	}

protected:
	TriMesh* mesh_;
	bool render_;
	SceneObjectId sceneId_;
	RenderMode renderMode_;
	bool rotationInvariant_;
};