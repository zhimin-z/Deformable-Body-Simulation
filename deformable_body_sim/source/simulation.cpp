// ---------------------------------------------------------------------------------//
// Copyright (c) 2013, Regents of the University of Pennsylvania                    //
// All rights reserved.                                                             //
//                                                                                  //
// Redistribution and use in source and binary forms, with or without               //
// modification, are permitted provided that the following conditions are met:      //
//     * Redistributions of source code must retain the above copyright             //
//       notice, this list of conditions and the following disclaimer.              //
//     * Redistributions in binary form must reproduce the above copyright          //
//       notice, this list of conditions and the following disclaimer in the        //
//       documentation and/or other materials provided with the distribution.       //
//     * Neither the name of the <organization> nor the                             //
//       names of its contributors may be used to endorse or promote products       //
//       derived from this software without specific prior written permission.      //
//                                                                                  //
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  //
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    //
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           //
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY               //
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       //
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     //
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      //
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       //
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    //
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     //
//                                                                                  //
// Contact Tiantian Liu (ltt1598@gmail.com) if you have any questions.              //
//----------------------------------------------------------------------------------//

#pragma warning( disable : 4996)

#include "simulation.h"
#include "timer_wrapper.h"

Simulation::Simulation()

{
}

Simulation::~Simulation()
{
    clearConstraints();
}

void Simulation::Reset()
{    
    m_external_force.resize(m_mesh->m_system_dimension);

    setupConstraints();

    m_selected_attachment_constraint = NULL;
}

void Simulation::Update()
{
	// LOOK: main simulation loop

    // update external force
    calculateExternalForce();

	// get a reference of current position and velocity
	VectorX& x = m_mesh->m_current_positions;
	VectorX& v = m_mesh->m_current_velocities;
	// substepping timestep for explicit/implicit euler
	ScalarType dt = m_h / m_iterations_per_frame;
	for (unsigned int it = 0; it != m_iterations_per_frame; ++it)
	{
		// update cloth
		switch (m_integration_method)
		{
		case INTEGRATION_EXPLICIT_EULER:
			integrateExplicitEuler(x, v, dt);
			break;
		case INTEGRATION_EXPLICIT_RK2:
			integrateExplicitRK2(x, v, dt);
			break;
		case INTEGRATION_EXPLICIT_RK4:
			integrateExplicitRK4(x, v, dt);
			break;
		case INTEGRATION_EXPLICIT_SYMPLECTIC:
			integrateExplicitSymplectic(x, v, dt);
			break;
		case INTEGRATION_IMPLICIT_EULER_BARAFF_WITKIN:
			integrateImplicitBW(x, v, dt);
			break;
		}
	}
	// pbd integration method
	if (m_integration_method == INTEGRATION_POSITION_BASED_DYNAMICS)
	{
		integratePBD(x, v, m_iterations_per_frame);
	}

	// collision detection and resolution
	std::vector<CollisionInstance> collisions;
	detectCollision(x, collisions);
	resolveCollision(x, v, collisions);

    // damping
    dampVelocity();
}

void Simulation::DrawConstraints(const VBO& vbos)
{
    for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
    {
        (*it)->Draw(vbos);
    }
}

ScalarType Simulation::TryToSelectAttachmentConstraint(const EigenVector3& p0, const EigenVector3& dir)
{
    ScalarType ray_point_dist;
    ScalarType min_dist = 100.0;
    AttachmentConstraint* best_candidate = NULL;

    bool current_state_on = false;
    for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
    {
        AttachmentConstraint* ac;
        if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
        {
            ray_point_dist = ((ac->GetFixedPoint()-p0).cross(dir)).norm();
            if (ray_point_dist < min_dist)
            {
                min_dist = ray_point_dist;
                best_candidate = ac;
            }
        }
    }
    // exit if no one fits
    if (min_dist > DEFAULT_SELECTION_RADIUS)
    {
        UnselectAttachmentConstraint();

        return -1;
    }
    else
    {
        SelectAtttachmentConstraint(best_candidate);
        EigenVector3 fixed_point_temp = m_mesh->m_current_positions.block_vector(m_selected_attachment_constraint->GetConstrainedVertexIndex());

        return (fixed_point_temp-p0).dot(dir); // this is m_cached_projection_plane_distance
    }
}

bool Simulation::TryToToggleAttachmentConstraint(const EigenVector3& p0, const EigenVector3& dir)
{
    EigenVector3 p1;

    ScalarType ray_point_dist;
    ScalarType min_dist = 100.0;
    unsigned int best_candidate = 0;
    // first pass: choose nearest point
    for (unsigned int i = 0; i != m_mesh->m_vertices_number; i++)
    {
        p1 = m_mesh->m_current_positions.block_vector(i);
        
        ray_point_dist = ((p1-p0).cross(dir)).norm();
        if (ray_point_dist < min_dist)
        {
            min_dist = ray_point_dist;
            best_candidate = i;
        }
    }
    for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
    {
        AttachmentConstraint* ac;
        if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
        {
            ray_point_dist = ((ac->GetFixedPoint()-p0).cross(dir)).norm();
            if (ray_point_dist < min_dist)
            {
                min_dist = ray_point_dist;
                best_candidate = ac->GetConstrainedVertexIndex();
            }
        }
    }
    // exit if no one fits
    if (min_dist > DEFAULT_SELECTION_RADIUS)
    {
        return false;
    }
    // second pass: toggle that point's fixed position constraint
    bool current_state_on = false;
    for (std::vector<Constraint*>::iterator c = m_constraints.begin(); c != m_constraints.end(); ++c)
    {
        AttachmentConstraint* ac;
        if (ac = dynamic_cast<AttachmentConstraint*>(*c)) // is attachment constraint
        {
            if (ac->GetConstrainedVertexIndex() == best_candidate)
            {
                current_state_on = true;
                m_constraints.erase(c);
                break;
            }
        }
    }
    if (!current_state_on)
    {
        AddAttachmentConstraint(best_candidate);
    }

    return true;
}

void Simulation::SelectAtttachmentConstraint(AttachmentConstraint* ac)
{
    m_selected_attachment_constraint = ac;
    m_selected_attachment_constraint->Select();
}

void Simulation::UnselectAttachmentConstraint()
{
    if (m_selected_attachment_constraint)
    {
        m_selected_attachment_constraint->UnSelect();
    }
    m_selected_attachment_constraint = NULL;
}

void Simulation::AddAttachmentConstraint(unsigned int vertex_index)
{
    AttachmentConstraint* ac = new AttachmentConstraint(&m_stiffness_attachment, &m_stiffness_attachment_pbd, vertex_index, m_mesh->m_current_positions.block_vector(vertex_index));
    m_constraints.push_back(ac);
}

void Simulation::MoveSelectedAttachmentConstraintTo(const EigenVector3& target)
{
    if (m_selected_attachment_constraint)
        m_selected_attachment_constraint->SetFixedPoint(target);
}

void Simulation::clearConstraints()
{
    for (unsigned int i = 0; i < m_constraints.size(); ++i)
    {
        delete m_constraints[i];
    }
    m_constraints.clear();
}

void Simulation::setupConstraints()
{
	// LOOK setup spring constraints
    clearConstraints();

    switch(m_mesh->m_mesh_type)
    {
    case MESH_TYPE_CLOTH:
        // procedurally generate constraints including to attachment constraints
        {
            // generating attachment constraints.
            AddAttachmentConstraint(0);
            AddAttachmentConstraint(m_mesh->m_dim[1]*(m_mesh->m_dim[0]-1));

			// TODO
            // generate stretch constraints. assign a stretch constraint for each edge.
			EigenVector3 p1, p2;
            for(std::vector<Edge>::iterator e = m_mesh->m_edge_list.begin(); e != m_mesh->m_edge_list.end(); ++e)
            {
                p1 = m_mesh->m_current_positions.block_vector(e->m_v1);
                p2 = m_mesh->m_current_positions.block_vector(e->m_v2);
                SpringConstraint* c = new SpringConstraint(&m_stiffness_stretch, &m_stiffness_stretch_pbd, e->m_v1, e->m_v2, (p1-p2).norm());
                m_constraints.push_back(c);
            }

			// TODO
            // generate bending constraints. naive solution using cross springs 

            for(int i = 0; i < m_mesh->m_dim[0]; ++i)
            {
                for(int j = 0; j < m_mesh->m_dim[1]; ++j)
                {
                    int p1_vect_index = m_mesh->m_dim[1] * i + j;
                    p1 = m_mesh->m_current_positions.block_vector(p1_vect_index);
					
					
                    if (i+2 < m_mesh->m_dim[0])
                    {
                        int p2_vect_index = (m_mesh->m_dim[1] * (i + 2)) + j;
                        p2 = m_mesh->m_current_positions.block_vector(p2_vect_index);
                        SpringConstraint* c = new SpringConstraint(&m_stiffness_bending, &m_stiffness_bending_pbd, p1_vect_index, p2_vect_index, (p1-p2).norm());
                        m_constraints.push_back(c);
                    }

					if (j+2 < m_mesh->m_dim[1])
                    {
                        int p2_vect_index = (m_mesh->m_dim[1] * i) + (j + 2);
                        p2 = m_mesh->m_current_positions.block_vector(p2_vect_index);
                        SpringConstraint* c = new SpringConstraint(&m_stiffness_bending, &m_stiffness_bending_pbd, p1_vect_index, p2_vect_index, (p1-p2).norm());
                        m_constraints.push_back(c);
                    }               
                }
            }
        
        }
        break;
    case MESH_TYPE_TET:
        {
            // generate stretch constraints. assign a stretch constraint for each edge.
            EigenVector3 p1, p2;
            for(std::vector<Edge>::iterator e = m_mesh->m_edge_list.begin(); e != m_mesh->m_edge_list.end(); ++e)
            {
                p1 = m_mesh->m_current_positions.block_vector(e->m_v1);
                p2 = m_mesh->m_current_positions.block_vector(e->m_v2);
                SpringConstraint *c = new SpringConstraint(&m_stiffness_stretch, &m_stiffness_stretch_pbd, e->m_v1, e->m_v2, (p1-p2).norm());
                m_constraints.push_back(c);
            }
        }
        break;
    }
}

void Simulation::dampVelocity()
{
    if (std::abs(m_damping_coefficient) < EPSILON)
        return;

    // post-processing damping
    EigenVector3 pos_mc(0.0, 0.0, 0.0), vel_mc(0.0, 0.0, 0.0);
    unsigned int i, size;
    ScalarType denominator(0.0), mass(0.0);
    size = m_mesh->m_vertices_number;
    for(i = 0; i < size; ++i)
    {
        mass = m_mesh->m_mass_matrix.coeff(i*3, i*3);

        pos_mc += mass * m_mesh->m_current_positions.block_vector(i);
        vel_mc += mass * m_mesh->m_current_velocities.block_vector(i);
        denominator += mass;
    }
    assert(denominator != 0.0);
    pos_mc /= denominator;
    vel_mc /= denominator;

    EigenVector3 angular_momentum(0.0, 0.0, 0.0), r(0.0, 0.0, 0.0);
    EigenMatrix3 inertia, r_mat;
    inertia.setZero(); r_mat.setZero();

    for(i = 0; i < size; ++i)
    {
        mass = m_mesh->m_mass_matrix.coeff(i*3, i*3);

        r = m_mesh->m_current_positions.block_vector(i) - pos_mc;
        angular_momentum += r.cross(mass * m_mesh->m_current_velocities.block_vector(i));

        //r_mat = EigenMatrix3(0.0,  r.z, -r.y,
        //                    -r.z, 0.0,  r.x,
        //                    r.y, -r.x, 0.0);

        r_mat.coeffRef(0, 1) = r[2];
        r_mat.coeffRef(0, 2) = -r[1];
        r_mat.coeffRef(1, 0) = -r[2];
        r_mat.coeffRef(1, 2) = r[0];
        r_mat.coeffRef(2, 0) = r[1];
        r_mat.coeffRef(2, 1) = -r[0];

        inertia += r_mat * r_mat.transpose() * mass;
    }
    EigenVector3 angular_vel = inertia.inverse() * angular_momentum;

    EigenVector3 delta_v(0.0, 0.0, 0.0);
    for(i = 0; i < size; ++i)
    {
        r = m_mesh->m_current_positions.block_vector(i) - pos_mc;
        delta_v = vel_mc + angular_vel.cross(r) - m_mesh->m_current_velocities.block_vector(i);     
        m_mesh->m_current_velocities.block_vector(i) += m_damping_coefficient * delta_v;
    }
}

void Simulation::calculateExternalForce()
{
	VectorX external_acceleration_field(m_mesh->m_system_dimension);
	external_acceleration_field.setZero();

    // LOOK: adding gravity
    for (unsigned int i = 0; i < m_mesh->m_vertices_number; ++i)
    {
        external_acceleration_field[3*i+1] += -m_gravity_constant;
    }

    m_external_force = m_mesh->m_mass_matrix * external_acceleration_field;
}

void Simulation::detectCollision(const VectorX& x, std::vector<CollisionInstance>& collisions)
{
    // Naive implementation of collision detection
    EigenVector3 normal;
    ScalarType dist;
	collisions.clear();

    for (unsigned int i = 0; i < m_mesh->m_vertices_number; ++i)
    {
        EigenVector3 xi = x.block_vector(i);

        if (m_scene->StaticIntersectionTest(xi, normal, dist)) 
        {
			// if collision
			collisions.push_back(CollisionInstance(i, normal, dist));
        }
    }
}
void Simulation::resolveCollision(VectorX&x, VectorX& v, std::vector<CollisionInstance>& collisions)
{
	for (std::vector<CollisionInstance>::iterator it = collisions.begin(); it!= collisions.end(); ++it)
	{
		// find the vertex id
		unsigned int id = it->point_index;
		// correct position
		x.block_vector(id) -= it->normal * it->dist;
		// correct velocity
		ScalarType vn = it->normal.dot(v.block_vector(id));
		v.block_vector(id) += -(1+m_restitution_coefficient)*vn*it->normal;
	}
}

void Simulation::integrateExplicitEuler(VectorX& x, VectorX& v, ScalarType dt)
{
	// TODO:
	// v_next = v_current + dt * a_current
	// x_next = x_current + dt * v_current

	VectorX current_force;
	computeForces(x, current_force);

	VectorX a = m_mesh->m_inv_mass_matrix * current_force;

	x = x + dt * v;
	v = v + dt * a;
	
	
}

void Simulation::integrateExplicitRK2(VectorX& x, VectorX& v, ScalarType dt)
{
	// TODO:
	// heun's method (RK2)

	ScalarType half_dt = dt / 2.0;
	VectorX temp_x, temp_v; // temp_x and temp_v are used to compute the forces in an intermediate state

	// step 1: compute k1
	VectorX f1;
	computeForces(x, f1);
	// get the slope for x and v separately 
	VectorX k1_x = v;
	VectorX k1_v = m_mesh->m_inv_mass_matrix * f1;
	// setup a temporary point for next iteration
	temp_x = x + half_dt * k1_x;
	temp_v = v + half_dt * k1_v;

	// step2: compute k2
	VectorX f2;
	computeForces(temp_x, f2);
	// get the slope for x and v separately 
	VectorX k2_x = temp_v;
	VectorX k2_v = m_mesh->m_inv_mass_matrix * f2;

	x = x + dt * k2_x;
	v = v + dt * k2_v;
}

void Simulation::integrateExplicitRK4(VectorX& x, VectorX& v, ScalarType dt)
{
	// LOOK: here's a sample work flow of RK4
	ScalarType half_dt = dt / 2.0;
	VectorX temp_x, temp_v; // temp_x and temp_v are used to compute the forces in an intermediate state

	// step 1: compute k1
	VectorX f1;
	computeForces(x, f1);
	// get the slope for x and v separately 
	VectorX k1_x = v;
	VectorX k1_v = m_mesh->m_inv_mass_matrix * f1;
	// setup a temporary point for next iteration
	temp_x = x + half_dt * k1_x;
	temp_v = v + half_dt * k1_v;

	// step2: compute k2
	VectorX f2;
	computeForces(temp_x, f2);
	// get the slope for x and v separately 
	VectorX k2_x = temp_v;
	VectorX k2_v = m_mesh->m_inv_mass_matrix * f2;
	// setup a temporary point for next iteration
	temp_x = x + half_dt * k2_x;
	temp_v = v + half_dt * k2_v;

	// step3: compute k3
	VectorX f3;
	computeForces(temp_x, f3);
	// get the slope for x and v separately 
	VectorX k3_x = temp_v;
	VectorX k3_v = m_mesh->m_inv_mass_matrix * f3;
	// setup a temporary point for next iteration
	temp_x = x + dt * k3_x;
	temp_v = v + dt * k3_v;

	// step4: compute k4
	VectorX f4;
	computeForces(temp_x, f4);
	// get the slope for x and v separately 
	VectorX k4_x = temp_v;
	VectorX k4_v = m_mesh->m_inv_mass_matrix * f4;

	// Put it all together
	x = x + 1/6.0 * dt * (k1_x + 2*k2_x + 2*k3_x + k4_x);
	v = v + 1/6.0 * dt * (k1_v + 2*k2_v + 2*k3_v + k4_v);
}

void Simulation::integrateExplicitSymplectic(VectorX& x, VectorX& v, ScalarType dt)
{
	// TODO:
	// v_next = v_current + dt * a_current
	// x_next = x_current + dt * v_next

	VectorX current_force;
	computeForces(x, current_force);

	VectorX a = m_mesh->m_inv_mass_matrix * current_force;

	v = v + dt * a;
	x = x + dt * v;
}

void Simulation::integrateImplicitBW(VectorX& x, VectorX& v, ScalarType dt)
{
	// TODO:
	// v_next = v_current + dt * a_next
	// x_next = x_current + dt * v_next
	// [M - dt^2 K] * v_next = M * v_current + dt * f(x_current);
	// A * v_next = b;

	VectorX current_force;
	computeForces(x, current_force);


	SparseMatrix stiffness_matrix;
    computeStiffnessMatrix(x, stiffness_matrix);

	ScalarType h_square = dt*dt;

	VectorX b = m_mesh->m_mass_matrix * v + dt * current_force;
    SparseMatrix A = m_mesh->m_mass_matrix - h_square * stiffness_matrix;

	Eigen::SimplicialLDLT<SparseMatrix, Eigen::Upper> ldltSolver;
    factorizeDirectSolverLDLT(A, ldltSolver, "BW method");
    VectorX v_next = ldltSolver.solve(b);

	v = v_next;
	x = x + dt * v;
}

void Simulation::integratePBD(VectorX& x, VectorX& v, unsigned int ns)
{
	// TODO:
	VectorX a = m_mesh->m_inv_mass_matrix * m_external_force;

	v = v + m_h * a;
	VectorX temp = x + m_h * v;

	bool change = true;

	for (unsigned int i = 0; i <= ns; ++i){
		if(change){
			for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
			{
				(*it)->PBDProject(temp, m_mesh->m_inv_mass_matrix, ns);
			}
			change = false;
		}else{
			for (std::vector<Constraint*>::reverse_iterator it = m_constraints.rbegin(); it != m_constraints.rend(); ++it)
			{
				(*it)->PBDProject(temp, m_mesh->m_inv_mass_matrix, ns);
			}
			change = true;
		}

	}

	v = (temp - x) / m_h;
	x = temp;
}

#pragma region implicit/explicit euler
void Simulation::computeForces(const VectorX& x, VectorX& force)
{
    VectorX gradient;

    gradient.resize(m_mesh->m_system_dimension);
    gradient.setZero();

    // springs
    for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
    {
        (*it)->EvaluateGradient(x, gradient);
    }

	// internal_force = - gradient of elastic energy
    force = -gradient;

	// external forces
	force += m_external_force;
}

void Simulation::computeStiffnessMatrix(const VectorX& x, SparseMatrix& stiffness_matrix)
{
	SparseMatrix hessian;
	hessian.resize(m_mesh->m_system_dimension, m_mesh->m_system_dimension);
	std::vector<SparseMatrixTriplet> h_triplets;
	h_triplets.clear();

	for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
	{
		(*it)->EvaluateHessian(x, h_triplets);
	}

	hessian.setFromTriplets(h_triplets.begin(), h_triplets.end());

	// stiffness_matrix = - hessian matrix
	stiffness_matrix = - hessian;
}

#pragma endregion

#pragma region utilities
void Simulation::factorizeDirectSolverLLT(const SparseMatrix& A, Eigen::SimplicialLLT<SparseMatrix, Eigen::Upper>& lltSolver, char* warning_msg)
{
    SparseMatrix A_prime = A;
    lltSolver.analyzePattern(A_prime);
    lltSolver.factorize(A_prime);
    ScalarType Regularization = 0.00001;
    bool success = true;
    while (lltSolver.info() != Eigen::Success)
    {
        Regularization *= 10;
        A_prime = A_prime + Regularization*m_mesh->m_identity_matrix;
        lltSolver.factorize(A_prime);
        success = false;
    }
    if (!success)
        std::cout << "Warning: " << warning_msg <<  " adding "<< Regularization <<" identites.(llt solver)" << std::endl;
}

void Simulation::factorizeDirectSolverLDLT(const SparseMatrix& A, Eigen::SimplicialLDLT<SparseMatrix, Eigen::Upper>& ldltSolver, char* warning_msg)
{
    SparseMatrix A_prime = A;
    ldltSolver.analyzePattern(A_prime);
    ldltSolver.factorize(A_prime);
    ScalarType Regularization = 0.00001;
    bool success = true;
    while (ldltSolver.info() != Eigen::Success)
    {
        Regularization *= 10;
        A_prime = A_prime + Regularization*m_mesh->m_identity_matrix;
        ldltSolver.factorize(A_prime);
        success = false;
    }
    if (!success)
        std::cout << "Warning: " << warning_msg <<  " adding "<< Regularization <<" identites.(ldlt solver)" << std::endl;
}

void Simulation::generateRandomVector(const unsigned int size, VectorX& x)
{
    x.resize(size);

    for (unsigned int i = 0; i < size; i++)
    {
        x(i) = ((ScalarType)(rand())/(ScalarType)(RAND_MAX+1)-0.5)*2;
    }
}

#pragma endregion
