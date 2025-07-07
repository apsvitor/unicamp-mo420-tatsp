#
# Este deve ser o unico arquivo que deve ser atualizado pelo aluno
#
# Abaixo voce encontra as rotinas vazias das seguintes funcoes:
#                 TriggerArcTSP_lb_lp(T)
#                 TriggerArcTSP_lb_rlxlag(T)
#                 TriggerArcTSP_lb_colgen(T) - Opcional
#                 TriggerArcTSP_ub_lp(T)
#                 TriggerArcTSP_ub_rlxlag(T)
#                 TriggerArcTSP_ub_colgen(T) - Opcional
#                 TriggerArcTSP_ilp(T)
#
using Dates
using Infiltrator
using JuMP
using Gurobi
using LinearAlgebra

# --------------------------------------------------------------
function TriggerArcTSP_lb_lp(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here



	T.time_lb_lp = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
function TriggerArcTSP_lb_rlxlag(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Dualize some of the ILP constraints.
	active_constraints = [
		true, 	# C1
		true,	# C2
		true, 	# C3
		true,	# C4
		true,	# C5
		false,	# C6
		true, #false,	# C7
		true, #false,	# C8
		false,	# C9
		true, #false	# C10
	]
	# Build the basic ILP model without the dualized constraints.
	model, model_vars = _build_ilp(T, active_constraints)

	# Initialize Lagrangean multipliers for the constraints.
	lambda_dict = Dict{String, Any}(
		"C1" => zeros(1),
		"C2" => zeros(T.NArcs),
		"C3" => zeros(T.NNodes),
		"C4" => zeros(T.NNodes),
		"C5" => zeros(T.NArcs),
		"C6" => zeros(T.NTriggers),
		"C7" => zeros(T.NTriggers),
		"C8" => zeros(T.NTriggers),
		"C9" => zeros(T.NTriggers),
		"C10" => zeros(T.NTriggers, T.NTriggers)
	)

	# Violation vectors for each dualized constraint set.
	violations_dict = Dict{String, Any}(
		"v1" => zeros(1),
		"v2" => zeros(T.NArcs),
		"v3" => zeros(T.NNodes),
		"v4" => zeros(T.NNodes),
		"v5" => zeros(T.NArcs),
		"v6" => zeros(T.NTriggers),
		"v7" => zeros(T.NTriggers),
		"v8" => zeros(T.NTriggers),
		"v9" => zeros(T.NTriggers),
		"v10" => zeros(T.NTriggers, T.NTriggers)
	)

	# Define bounds
	max_arc_cost = [a.cost for a in T.Arc]
	for r in T.Trigger
		target_arc = r.target_arc_id
		max_arc_cost[target_arc] = max(max_arc_cost[target_arc], r.cost)
	end

	# This is not necessarily a feasible cost but it guarantees that it is the
	# maximum tour cost if the most expensive arcs were traversed and their activations
	# were the most expensive ones possible.
	UB = sum(sort(max_arc_cost, rev=true)[1:T.NNodes])
	LB = -Inf

	# Stores the best Lagrangean multipliers found so far
	# Needs to be a deepcopy because lambda_dict will be updated in each iteration
	best_lambda = deepcopy(lambda_dict)
	best_lambda = 
	z_k_values = []

	set_silent(model)
	MAX_ITERS = 20
	max_time_per_iteration = floor(T.maxtime_lb_rlxlag / MAX_ITERS)
	set_attribute(model, "TimeLimit", max_time_per_iteration)
	# Subgradient method
	for k in 1:MAX_ITERS
		# Update the model with the current penalization values
		model = _update_lag_obj_function(
			T,
			active_constraints,
			model,
			model_vars,
			lambda_dict
		)

		optimize!(model)

		# Check the optimization status
		_optimization_status_check(model)

		# Get the current objective value and updates LB if it is better
		z_k = objective_value(model)
		if z_k > LB
			LB = z_k
			best_lambda = deepcopy(lambda_dict)
			println("Improved LB at iteration[$(k)]: $LB")
		end
		# Storing for debugging purposes (might delete later)
		push!(z_k_values, z_k)

		# Get the new variable values
		current_vars = Dict{String, Any}(
			"x" => JuMP.value.(model_vars["x"]),
			"y" => JuMP.value.(model_vars["y"]),
			"y_r" => JuMP.value.(model_vars["y_r"]),
			"y_hat" => JuMP.value.(model_vars["y_hat"]),
			"u" => JuMP.value.(model_vars["u"])
		)

		violations_dict = _compute_lag_violations(
			T,
			active_constraints,
			current_vars,
			violations_dict
		)
		# println("C6: $([v6 for v6 in violations_dict["v6"] if v6 > 0])")
		# println("C9: $([v9 for v9 in violations_dict["v9"] if v9 > 0])")

		# I'll add an extra stopping criterion based on violations
		v_norm = norm(vcat(
			violations_dict["v1"],
			violations_dict["v2"],
			violations_dict["v3"],
			violations_dict["v4"],
			violations_dict["v5"],
			violations_dict["v6"],
			violations_dict["v7"],
			violations_dict["v8"],
			violations_dict["v9"],
			vec(violations_dict["v10"])
			), 2)
		if v_norm <= 1e-3
			println("Stopping criteria met at iteration $(k): violations are small enough.")
			break
		end
		alpha_k = 0
		beta = 0.7
		# Polyak's currently not doing great (probably my fault T.T)
		if v_norm > 0
			alpha_k = min(5, beta * (UB - z_k) / v_norm^2)
		end
		alpha_k = 2.0/k
		println("alpha_$(k) = $(alpha_k) | v_norm2 = $(v_norm^2)")

		lambda_dict = _update_lag_multipliers(
			T,
			lambda_dict,
			violations_dict,
			alpha_k
		)

		if ceil(time() - StartingTime) >= T.maxtime_lb_rlxlag
			println("Stopping criteria met at iteration $(k): time limit reached.")
			break
		end

	end
	println("Lagrangean relaxation LB: $(LB) | UB: $(UB)")
	println(z_k_values)
	_optimization_statistics(model)
	_print_path(model_vars)
	T.lb_rlxlag = LB
	T.time_lb_rlxlag = ceil(time() - StartingTime) # To compute the time used by this routine.
end

function _update_lag_obj_function(
	T::TriggerArcTSP,
	active_constraints::Vector{Bool},
	model::Model,
	model_vars::Dict{String, Any},
	lambda_dict::Dict{String, Any})
	"""
	Updates the objective function of the model to include the Lagrangean multipliers.
	Args:
		T (TriggerArcTSP): The TriggerArcTSP instance containing the problem data.
		model (Model): The JuMP model to be updated.
		model_vars (Dict{String, Any}): A dictionary containing the model variables.
			It should contain keys "x", "y", "y_r", "y_hat", and "u".
		lambda_dict (Dict{String, Vector{Float64}}): A dictionary containing the
			Lagrangean multipliers for the constraints C6 to C10.
	"""
	x = model_vars["x"]
	y = model_vars["y"]
	y_r = model_vars["y_r"]
	y_hat = model_vars["y_hat"]
	u = model_vars["u"]
	M = T.NNodes

	base_obj_terms = (
		sum(T.Trigger[r].cost * y_r[r] for r in 1:T.NTriggers)
		+ sum(T.Arc[a].cost * y[a] for a in 1:T.NArcs)
	)
	C1_term = if !active_constraints[1]
		lambda_dict["C1"] * (sum(x[a] for a in 1:T.NArcs) - T.NNodes)
	else
		0
	end
	C2_term = if !active_constraints[2]
		sum(lambda_dict["C2"][a] * (
			u[T.Arc[a].u]
			+ 1
			- u[T.Arc[a].v]
			- M * (1 - x[a])
			) for a in 1:T.NArcs if T.Arc[a].v != 1
		)
	else
		0
	end
	C3_term = if !active_constraints[3]
		sum(lambda_dict["C3"][v] * (
			sum(x[a] for a in T.NArcs if T.Arc[a].v == v; init=0) - 1
			) for v in 1:T.NNodes
		)
	else
		0
	end
	C4_term = if !active_constraints[4]
		sum(lambda_dict["C4"][v] * (
			sum(x[a] for a in T.NArcs if T.Arc[a].u == v; init=0) - 1
			) for v in 1:T.NNodes
		)
	else
		0
	end
	C5_term = if !active_constraints[5]
		sum(lambda_dict["C5"][a] * (
			y[a] + sum(y_r[r] for r in findall(t_a -> t_a.target_arc_id == a, T.Trigger); init=0) - x[a]
			) for a in 1:T.NArcs
		)
	else
		0
	end
	C6_term = if !active_constraints[6]
		sum(
			lambda_dict["C6"][r] * (
				y_r[r] - x[T.Trigger[r].trigger_arc_id]
			) for r in 1:T.NTriggers
		)
	else
		0
	end
	C7_term = if !active_constraints[7]
		sum(lambda_dict["C7"][r] * (
			u[T.Arc[T.Trigger[r].trigger_arc_id].u]
			+ 1
			- u[T.Arc[T.Trigger[r].target_arc_id].u]
			- M*(1 - y_r[r])
			) for r in 1:T.NTriggers
		)
	else
		0
	end
	C8_term = if !active_constraints[8]
		sum(lambda_dict["C8"][r] * (
			u[T.Arc[T.Trigger[r].target_arc_id].u]
			+ 1
			- u[T.Arc[T.Trigger[r].trigger_arc_id].u]
			- M*(1 - y_hat[r])
			) for r in 1:T.NTriggers
		)
	else
		0
	end
	C9_term = if !active_constraints[9]
		sum(lambda_dict["C9"][r] * (
			x[T.Trigger[r].trigger_arc_id]
			- (1 - x[T.Trigger[r].target_arc_id])
			- (1 - y[T.Trigger[r].target_arc_id])
			- y_hat[r]
			) for r in 1:T.NTriggers
		)
	else
		0
	end
	C10_term = if !active_constraints[10]
		sum(
			lambda_dict["C10"][r1, r2] * (
				u[T.Arc[T.Trigger[r2].trigger_arc_id].u]
				- M*y_hat[r2]
				- u[T.Arc[T.Trigger[r1].trigger_arc_id].u]
				- M*(2 - y_r[r1] - x[T.Trigger[r2].trigger_arc_id])
				+ 1
				) for r1 in 1:T.NTriggers for r2 in r1+1:T.NTriggers 
				if T.Trigger[r1].target_arc_id == T.Trigger[r2].target_arc_id
		)
	else
		0
	end

	@objective(
		model,
		Min,
		base_obj_terms
		+ C1_term
		+ C2_term
		+ C3_term
		+ C4_term
		+ C5_term
		+ C6_term
		+ C7_term
		+ C8_term
		+ C9_term
		+ C10_term
	)
	return model
end

function _compute_lag_violations(
	T::TriggerArcTSP,
	active_constraints::Vector{Bool},
	model_vars::Dict{String, Any},
	violations_dict::Dict{String, Any},
	)
	"""Computes the violations for the Lagrangean multipliers.

	Args:
		T (TriggerArcTSP): The TriggerArcTSP instance containing the problem data.
		model_vars (Dict{String, Any}): A dictionary containing the model variables.
			It should contain keys "x", "y", "y_r", "y_hat", and "u".
		lambda_dict (Dict{String, Vector{Float64}}): A dictionary containing the
			Lagrangean multipliers for the constraints C6 to C10.
		violations_dict (Dict{String, Vector{Float64}}): A dictionary to store the
			violations for each constraint.
	Returns:
		violations_dict: dictionary of violations for each constraint.
	"""
	x = model_vars["x"]
	y = model_vars["y"]
	y_r = model_vars["y_r"]
	y_hat = model_vars["y_hat"]
	u = model_vars["u"]
	M = T.NNodes

	# Apply a scaling method to avoid big violations related to the difference
	# between BigM and the original costs from the objetive function.
	c_max = maximum([a.cost for a in T.Arc])
	for r in T.Trigger
		c_max = max(c_max, r.cost)
	end
	scale = max(c_max, M, 1.0)
	# Deactivate scaling
	scale = 1.0

	# C1
	if active_constraints[1] == false
		violations_dict["v1"][0] = sum(x[a] for a in 1:T.NArcs) - T.NNodes
	end

	# C2
	if active_constraints[2] == false
		for a in 1:T.NArcs
			i = T.Arc[a].u
			j = T.Arc[a].v
			if j != 1
				violations_dict["v2"][a] = max(0, u[i] + 1 - u[j] - M * (1 - x[a]))
			end
		end
	end

	# C3
	if active_constraints[3] == false
		for v in 1:T.NNodes
			violations_dict["v3"][v] = sum(x[a] for a in 1:T.NArcs if T.Arc[a].v == v) - 1
		end
	end

	# C4
	if active_constraints[4] == false
		for v in 1:T.NNodes
			violations_dict["v4"][v] = sum(x[a] for a in 1:T.NArcs if T.Arc[a].u == v) - 1
		end
	end

	# C5
	if active_constraints[5] == false
		for a in 1:T.NArcs
			R_a = findall(triggers_a -> triggers_a.target_arc_id == a, T.Trigger)
			violations_dict["v5"][a] = y[a] + sum(y_r[r] for r in R_a) - x[a]
		end
	end

	# C6
	if active_constraints[6] == false
		for r in 1:T.NTriggers
			trigger = T.Trigger[r].trigger_arc_id
			violations_dict["v6"][r] = max(0, y_r[r] - x[trigger])
		end
	end

	# C7
	if active_constraints[7] == false
		for r in 1:T.NTriggers
			i = T.Arc[T.Trigger[r].trigger_arc_id].u
			h = T.Arc[T.Trigger[r].target_arc_id].u
			violations_dict["v7"][r] = max(0, u[i] + 1 - u[h] - M*(1 - y_r[r])) / scale
		end
	end

    # C8
	if active_constraints[8] == false
		for r in 1:T.NTriggers
			i = T.Arc[T.Trigger[r].trigger_arc_id].u
			h = T.Arc[T.Trigger[r].target_arc_id].u
			violations_dict["v8"][r] = max(0, u[h] + 1 - u[i] - M*(1 - y_hat[r])) / scale
		end
	end

    # C9
	if active_constraints[9] == false
		for r in 1:T.NTriggers
			trigger = T.Trigger[r].trigger_arc_id
			target  = T.Trigger[r].target_arc_id
			violations_dict["v9"][r] = max(0,
				x[trigger]
				- (1 - x[target])
				- (1 - y[target])
				- y_hat[r]
			)
		end
	end

    # C10
	if active_constraints[10] == false
		for r1 in 1:T.NTriggers
			for r2 in r1+1:T.NTriggers
				if T.Trigger[r1].target_arc_id == T.Trigger[r2].target_arc_id
					i_r1   = T.Arc[T.Trigger[r1].trigger_arc_id].u
					i_r2   = T.Arc[T.Trigger[r2].trigger_arc_id].u
					a_hat2 = T.Trigger[r2].trigger_arc_id
					violations_dict["v10"][r1,r2] = max(0,
						u[i_r2]
						- M*y_hat[r2]
						- u[i_r1]
						- M*(2 - y_r[r1] - x[a_hat2])
						+ 1
					) / scale
				end
			end
		end
	end

	return violations_dict
end

function _update_lag_multipliers(
	T::TriggerArcTSP,
	lambda_dict::Dict{String, Any},
	violations_dict::Dict{String, Any},
	alpha_k::Float64)
	"""Updates the Lagrangean multipliers based on the violations.

	Args:
		T(TriggerArcTSP): The TriggerArcTSP instance containing the problem data.
		lambda_dict (Dict{String, Vector{Float64}}): A dictionary containing the
			Lagrangean multipliers for the constraints C6 to C10.
		violations_dict (Dict{String, Vector{Float64}}): A dictionary containing the
			violations for each constraint.
		alpha_k (Float64): The step size for the update of the multipliers.
			It is computed based on the current violations and the Lagrangean objective.
	Returns:
		lambda_dict: Updated dictionary of Lagrangean multipliers.
	"""
	for v in 1:T.NNodes
		lambda_dict["C3"][v] = lambda_dict["C3"][v] + alpha_k * violations_dict["v3"][v]
		lambda_dict["C4"][v] = lambda_dict["C4"][v] + alpha_k * violations_dict["v4"][v]
	end

	for a in 1:T.NArcs
		lambda_dict["C2"][a] = lambda_dict["C2"][a] + alpha_k * violations_dict["v2"][a]
		lambda_dict["C5"][a] = lambda_dict["C5"][a] + alpha_k * violations_dict["v5"][a]
	end

	# I could use r1 loop at C10 to compute these guys but is it worth it?
	for r in 1:T.NTriggers
		lambda_dict["C6"][r] = lambda_dict["C6"][r] + alpha_k * violations_dict["v6"][r]
		lambda_dict["C7"][r] = lambda_dict["C7"][r] + alpha_k * violations_dict["v7"][r]
		lambda_dict["C8"][r] = lambda_dict["C8"][r] + alpha_k * violations_dict["v8"][r]
		lambda_dict["C9"][r] = lambda_dict["C9"][r] + alpha_k * violations_dict["v9"][r]
	end

	# maybe saving some time by never updating multipliers for unnecessary pairs
	for r1 in 1:T.NTriggers
		for r2 in r1+1:T.NTriggers
			if T.Trigger[r1].target_arc_id == T.Trigger[r2].target_arc_id
				lambda_dict["C10"][r1, r2] = lambda_dict["C10"][r1, r2] + alpha_k * violations_dict["v10"][r1, r2]
			end
		end
	end

	return lambda_dict
end

# --------------------------------------------------------------
function TriggerArcTSP_lb_colgen(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here




	T.time_lb_colgen = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
function TriggerArcTSP_ub_lp(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here




	T.time_ub_lp = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
function TriggerArcTSP_ub_rlxlag(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here




	T.time_ub_rlxlag = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
function TriggerArcTSP_ub_colgen(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here




	T.time_ub_colgen = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
function TriggerArcTSP_ilp(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	model, model_vars = _build_ilp(T)
	set_attribute(model, "TimeLimit", T.maxtime_ilp)


	optimize!(model)

	_optimization_statistics(model)
	_optimization_status_check(model)
	_print_path(model_vars)

	# To compute the time used by this routine.
	T.time_ilp = ceil(time() - StartingTime)
end
# --------------------------------------------------------------


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Auxiliary routines to simplify and avoid redundances
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function _build_ilp(
	T::TriggerArcTSP,
	constraint_selection::Vector{Bool}=fill(true, 10))
	"""Builds the ILP model for the Trigger Arc TSP problem.

	Args:
		T (TriggerArcTSP): The TriggerArcTSP instance containing the problem data.
		constraint_selection (Vector{Bool}): A vector indicating which constraints to
			include in the model.
			Each element corresponds to a constraint (from C1 to C10, C0 is mandatory),
			where `true` means the constraint is included, and `false` means it is not.
	Returns:
		Model: The JuMP model representing the ILP.
	"""
	model = Model(Gurobi.Optimizer)
	
	# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# # 							Model Variables
	# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# [x_a] -> x_ij: 1 if the arc (i, j) is in the solution, 0 otherwise.
	x = @variable(model, x[i=1:T.NArcs], Bin)

	# [y_a]: 1 if the arc a = (i, j) is in the solution and has no active relationship,
	# 0 otherwise.
	# Implication: if y_a = 1, then the cost of the arc will be the base cost: c(a).
	y = @variable(model, y[i=1:length(T.Arc)], Bin)

	for a in 1:T.NArcs
		set_name(x[a], "x_"*string(a))
		set_name(y[a], "y_"*string(a))
	end

	# [y_r]: 1 if the relationship r is active, 0 otherwise.
	y_r = @variable(model, y_r[1:T.NTriggers], Bin)
	
	# [y^_r] -> y_hat_r: 1 if the relationship r does not occur either because the
	# trigger arc occurred after the target arc or if the target arc is not in the
	# solution.
	y_hat = @variable(model, y_hat[1:T.NTriggers], Bin)
	
	for r in 1:T.NTriggers
		set_name(y_r[r], "y_r"*string(r))
		set_name(y_hat[r], "y^_r"*string(r))
	end

	# [u_i]: auxiliary variable indicating the visitation order of node i.
	u = @variable(model, u[i=1:T.NNodes], lower_bound=0, upper_bound=T.NNodes, Int)

	# Build a dictionary to store the variables for easy access later on.
	model_vars = Dict{String, Any}(
		"x" => x,
		"y" => y,
		"y_r" => y_r,
		"y_hat" => y_hat,
		"u" => u
	)

	# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# #								Objective Function
	# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# Minimize the total cost of the tour considering:
	# - The base cost of the arcs when they are not affected by any trigger.
	# - The modified cost of the arcs when they are affected by a trigger.
	@objective(
		model,
		Min,
		sum(T.Trigger[r].cost * y_r[r] for r in 1:T.NTriggers)
		+ sum(T.Arc[a].cost * y[a] for a in 1:T.NArcs)
	)

	# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# #								Model Constraints
	# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# C0: Starting node constraint (mandatory).
	__add_constraint_C0(model, u)

	# C1: Hamiltonian Cycle Constraint
	# DISCLAIMER: This constraint is redundant and will remain deactivated.
	# C3 & C4 implicitly ensure that the solution is a Hamiltonian cycle.
	if constraint_selection[1]
		__add_constraint_C1(model, x, T)
	end

	# C2: Visitation order of nodes.
	if constraint_selection[2]
		__add_constraint_C2(model, u, x, T)
	end

	# C3 & C4: In-degree and Out-degree Constraints, respectively.
	# Each node must have exactly one incoming and one outgoing arc.
	if constraint_selection[3]
		__add_constraint_C3(model, x, T)
	end
	if constraint_selection[4]
		__add_constraint_C4(model, x, T)
	end

	# C5: Ensures that the costs of the arcs are either the base cost or the modified
	# cost, and only if the arc is in the tour.
	if constraint_selection[5]
		__add_constraint_C5(model, y, y_r, x, T)
	end

	# C6: For each relation r = ((i,j),(h,k)), the corresponding variable y_r be active
	# if the tour traverses the trigger arc (i,j), i.e., if x_ij = 1
	if constraint_selection[6]
		__add_constraint_C6(model, y_r, x, T)
	end

	# C7: For each relation r = ((i,j),(h,k)) ensure that if the arc (i,j) does not
	# precede the arc (h,k) in the solution, then the variable y_r = 0.
	if constraint_selection[7]
		__add_constraint_C7(model, u, y_r, T)
	end

	# C8: If arc (h,k) does not precede arc (i,j) in the solution, then variable
	# y^r = 0, implying that variables y^r may only assume the value 1 if arc (h,k)
	# precedes arc (i,j) in the solution.
	if constraint_selection[8]
		__add_constraint_C8(model, u, y_hat, T)
	end

	# C9: For each relation r = ((i,j),(h,k)), the trigger arc (i,j) can only be used
	# if at least one of the following three conditions is satisfied:
	#	(1) arc (h,k) is not used;
	#	(2) variable y_hk associated with the base cost of arc (h,k) is not activated;
	#	(3) if both arcs (i,j) and (h,k) are used in the solution, the arc (i,j) is
	#		traversed after the arc (h,k).
	if constraint_selection[9]
		__add_constraint_C9(model, x, y, y_hat, T)
	end

	# C10: For each pair of relations r1 = ((i,j),(h,k)), r2 = ((^i, ^j), (h, k))
	# differing only by the trigger arc requires that the active relation is the one
	# whose trigger arc is traversed last but before the arc (h,k) in the tour.
	if constraint_selection[10]
		__add_constraint_C10(model, u, y_hat, y_r, x, T)
	end

	return model, model_vars
end

function __add_constraint_C0(
	model::Model,
	u::Vector{VariableRef})
	# C0: Starting node constraint.
	@constraint(
		model,
		u[1] == 1,
		base_name="C0"
	)
end

function __add_constraint_C1(
	model::Model,
	x::Vector{VariableRef},
	T::TriggerArcTSP)
	@constraint(
		model,
		sum(x[a] for a in 1:T.NArcs) == T.NNodes,
		base_name="C1"
	)
end

function __add_constraint_C2(
	model::Model,
	u::Vector{VariableRef},
	x::Vector{VariableRef},
	T::TriggerArcTSP)
	BIG_M_C2 = T.NNodes
	for a in 1:T.NArcs
		i = T.Arc[a].u
		j = T.Arc[a].v
		if j != 1
			@constraint(
				model,
				u[i] + 1 <= u[j] + BIG_M_C2 * (1 - x[a]),
				base_name="C2"
			)
		end
	end
end

function __add_constraint_C3(
	model::Model,
	x::Vector{VariableRef},
	T::TriggerArcTSP)
	for v in 1:T.NNodes
		@constraint(
			model,
			sum(x[a] for a in 1:T.NArcs if T.Arc[a].v == v) == 1,
			base_name="C3"
		)
	end
end

function __add_constraint_C4(
	model::Model,
	x::Vector{VariableRef},
	T::TriggerArcTSP)
	for v in 1:T.NNodes
		@constraint(
			model,
			sum(x[a] for a in 1:T.NArcs if T.Arc[a].u == v) == 1,
			base_name="C4"
		)
	end
end

function __add_constraint_C5(
	model::Model,
	y::Vector{VariableRef},
	y_r::Vector{VariableRef},
	x::Vector{VariableRef},
	T::TriggerArcTSP)
	for a in 1:T.NArcs
		# Set of relationships associated with arc a.
		R_a = findall(triggers_a -> triggers_a.target_arc_id == a, T.Trigger)

		@constraint(
			model,
			y[a] + sum(y_r[r] for r in R_a) == x[a],
			base_name="C5"
		)
	end
end

function __add_constraint_C6(
	model::Model,
	y_r::Vector{VariableRef},
	x::Vector{VariableRef},
	T::TriggerArcTSP)
	for r in 1:T.NTriggers
		trigger_arc = T.Trigger[r].trigger_arc_id
		@constraint(
			model,
			y_r[r] <= x[trigger_arc],
			base_name="C6"
		)
	end
end

function __add_constraint_C7(
	model::Model,
	u::Vector{VariableRef},
	y_r::Vector{VariableRef},
	T::TriggerArcTSP)
	BIG_M_C7 = T.NNodes
	for r in 1:T.NTriggers
		i = T.Arc[T.Trigger[r].trigger_arc_id].u
		h = T.Arc[T.Trigger[r].target_arc_id].u
		@constraint(
			model,
			u[i] + 1 <= u[h] +BIG_M_C7 * (1 - y_r[r]),
			base_name="C7"
		)
	end
end

function __add_constraint_C8(
	model::Model,
	u::Vector{VariableRef},
	y_hat::Vector{VariableRef},
	T::TriggerArcTSP)
	BIG_M_C8 = T.NNodes
	for r in 1: T.NTriggers
		h = T.Arc[T.Trigger[r].target_arc_id].u
		i = T.Arc[T.Trigger[r].trigger_arc_id].u
		@constraint(
			model,
			u[h] + 1 <= u[i] + BIG_M_C8 * (1 - y_hat[r]),
			base_name="C8"
		)
	end
end

function __add_constraint_C9(
	model::Model,
	x::Vector{VariableRef},
	y::Vector{VariableRef},
	y_hat::Vector{VariableRef},
	T::TriggerArcTSP)
	for r in 1:T.NTriggers
		trigger_arc = T.Trigger[r].trigger_arc_id
		target_arc = T.Trigger[r].target_arc_id
		@constraint(
			model,
			x[trigger_arc] <= (1 - x[target_arc]) + (1 - y[target_arc]) + y_hat[r],
			base_name="C9"
		)
	end
end

function __add_constraint_C10(
	model::Model,
	u::Vector{VariableRef},
	y_hat::Vector{VariableRef},
	y_r::Vector{VariableRef},
	x::Vector{VariableRef},
	T::TriggerArcTSP)
	BIG_M_C10 = T.NNodes
	for r1 in 1:T.NTriggers
		for r2 in r1+1:T.NTriggers
			if T.Trigger[r1].target_arc_id == T.Trigger[r2].target_arc_id
				a_hat_r2 = T.Trigger[r2].trigger_arc_id
				i_hat_r2 = T.Arc[a_hat_r2].u
				i_r1 = T.Arc[T.Trigger[r1].trigger_arc_id].u
				@constraint(
					model,
					u[i_hat_r2] - BIG_M_C10 * y_hat[r2] <= u[i_r1] + BIG_M_C10 * (2 - y_r[r1] - x[a_hat_r2]) - 1,
					base_name="C10"
				)
			end
		end
	end
end


function _optimization_status_check(model::Model)
	"""Checks the optimization status and handles problems.

	When a model is not optimal, computes conflicts to identify possible issues.
	Args:
		model (Model): The JuMP model to check.
	"""
	if termination_status(model) != MOI.OPTIMAL
		if termination_status(model) == MOI.TIME_LIMIT
			println("It was not possible to obtain optimal solution, due to time limit.")
			return MOI.TIME_LIMIT
		else
			compute_conflict!(model)
			if get_attribute(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
				iis_model, _ = copy_conflict(model)
				print(iis_model)
			end
			println(termination_status(model))
			error("It was not possible to obtain optimal solution.")
		end
	else
		println("Model solved successfully with status: ", termination_status(model))
		return MOI.OPTIMAL
	end
end

function _optimization_statistics(model::Model)
	"""Prints the optimization statistics of the model.

	Args:
		model (Model): The JuMP model to print statistics for.
	"""
	println("Optimization Statistics:")
	println("  Objective Value: ", objective_value(model))
	println("  Solve Time: ", solve_time(model), " seconds")
	println("  Number of Variables: ", num_variables(model))
	println("  Number of Constraints: ", num_constraints(model, count_variable_in_set_constraints=false))
	println("  Termination Status: ", termination_status(model))
end

function _print_path(model_vars::Dict{String, Any})
	"""Prints the path of the solution.

	Args:
		model (Model): The JuMP model containing the solution.
		T (TriggerArcTSP): The TriggerArcTSP instance containing the problem data.
	"""
	println("Solution Path:")
	tour = []
	u_values = model_vars["u"]
	# @infiltrate
	for i in 1:length(u_values)
		push!(tour, (i-1, Int(value(u_values[i])) ) )
	end
	sort!(tour, by=elem->elem[2])
	for i in tour
		print(i[1], ",")
	end
end