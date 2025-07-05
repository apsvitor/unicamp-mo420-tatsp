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

# --------------------------------------------------------------
function TriggerArcTSP_lb_lp(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.

	# Fill your routine here



	T.time_lb_lp = ceil(time() - StartingTime) # To compute the time used by this routine.
end
# --------------------------------------------------------------
function TriggerArcTSP_lb_rlxlag(T::TriggerArcTSP)
	StartingTime = time()  # To compute the time used by this routine.
	active_constraints = [
		true, 	# C1
		true,	# C2
		true, 	# C3
		true,	# C4
		true,	# C5
		false,	# C6
		false,	# C7
		false,	# C8
		false,	# C9
		false	# C10
	]
	model = _build_ilp(T, active_constraints)

	# First attempt: dualize some of the ILP constraints.
	# C6 to C10 seem to be the most "complicated" constraints.








	T.time_lb_rlxlag = ceil(time() - StartingTime) # To compute the time used by this routine.
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

	model = _build_ilp(T)

	optimize!(model)

	_optimization_statistics(model)
	_optimization_status_check(model)
	_print_path(model, T)

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
	# if constraint_selection[1]
		# __add_constraint_C1(model, x, T)
	# end

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

	return model
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
			error("It was not possible to obtain optimal solution, due to time limit.")
		else
			compute_conflict!(model)
			if get_attribute(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
				iis_model, _ = copy_conflict(model)
				print(iis_model)
			end
			println(termination_status(model))
			error("It was not possible to obtain optimal solution.")
		end
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

function _print_path(model::Model ,T::TriggerArcTSP)
	"""Prints the path of the solution.

	Args:
		model (Model): The JuMP model containing the solution.
		T (TriggerArcTSP): The TriggerArcTSP instance containing the problem data.
	"""
	println("Solution Path:")
	u_values = value.(model[:u])
	for v in u_values
		print("[", Int(v), "] -> ")
	end
	println("[", Int(u_values[1]), "]")
end