('est "Empirial System Theory (EST)"						;;;; Reference Tree

  ('System									;;; System Argument

    ('function-handle  "Vector field")
    ('function-handle  "Output functional")
    ('function-handle  "Adjoint vector field")
    ('positive-integer "Number of inputs")
    ('positive-integer "Number of States")
    ('positive-integer "Number of Outputs")
    ('positive-real    "Time step width")
    ('positive-real    "Time horizon"))

  ('Task									;;; Task Argument

    ('matrix_equation "Matrix equation solver"					;; Matrix equation task

      ('lyapunov  "Lyapunov Equation")
      ('sylvester "Sylvester Equation"))

    ('singular_values "Singular values of system gramian"			;; Gramian singular value task

      ('controllability "Controllability Gramian")
      ('observability   "Observability Gramian")
      ('minimality      "Cross Gramian"))

    ('model_reduction "Model reduction of state space"				;; Model reduction task

      ('poor_man            "Poor man aka POD")
      ('dominant_subspaces  "Dominant Subspaces aka DSPMR")
      ('approx_balancing    "Approximate balancing aka Modified POD")
      ('balanced_pod        "Balanced POD")
      ('balanced_truncation "Balanced truncation")
      ('dmd_galerkin        "Dynamic mode decomposition Galerkin")

        (('controllability "Controllability gramian only")
         ('observability   "Controllability and observability gramian")
         ('minimality      "Cross Gramian")))

    ('parameter_reduction "Parameter reduction"				;; Parameter reduction task

      ('observability "Identifiability gramian")
      ('minimality    "Cross-identifiability gramian"))

    ('combined_reduction "Combined state and parameter reduction"		;; Combined reduction task

      ('poor_man            "Poor man")
      ('dominant_subspaces  "Dominant subspaces")
      ('approx_balancing    "Approximate balancing")
      ('balanced_pod        "Balanced POD")
      ('balanced_truncation "Balanced truncation")

        (('observability "Controllability, observability and identifiability gramian")
         ('minimality    "Cross and cross-identifiability gramian")))

    ('decentralized_control "Decentralized Control"				;; Decentralized control task

      ('relative_gain_array  "Relative gain array")
      ('io_coherence         "Input-output coherence")
      ('io_pairing           "Input-output pairing")
      ('hardy_2              "Hardy-2 norm based")
      ('hardy_inf            "Hardy-Infinity norm based")
      ('participation_matrix "Participation matrix")
      ('hankel_interaction   "Hankel interaction array")
      ('rms_hsv              "Hilbert-Schmidt-Hankel norm based"))

    ('state_sensitivity "State sensitivity"					;; State sensitivity task

      ('controllability "Controllability-based")
      ('observability   "Observability-based")
      ('minimality      "Minimality-based"))

    ('parameter_sensitivity "Parameter sensitivity"				;; Parameter sensitivity task

      ('controllability "Controllability-based")
      ('observability   "Observability-based")
      ('minimality      "Minimality-based"))

    ('parameter_identifiability "Parameter identifiability"			;; Parameter identification task

      ('observability "Observability-based")
      ('minimality    "Minimality-based"))

    ('uncertainty_quantification "Uncertainty Quantification"			;; Uncertainty quantification task

      ('controllability "Controllability-based")
      ('observability   "Observability-based"))

    ('nonlinearity_quantification "Nonlinearity Quantification"		;; Nonlinearity quantification task

      ('controllability "Controllability-based")
      ('observability   "Observability-based")
      ('minimality      "Minimality-based")
      ('correlation     "Correlation-based"))

    ('gramian_index "Gramian index"						;; Gramian index task

      ('sigma_min             "Minimal singular value")			; -∞ Generalized mean
      ('harmonic_mean         "Harmonic mean of singular values")		; -1 Generalized mean
      ('geometric_mean        "Geometric mean of singular values")		; 0 Generalized mean
      ('energy_fraction       "Relative nuclear norm aka trace norm")		; 1 Generalized mean
      ('operator_norm         "Relative operator norm aka Hilbert-Schmidt norm"); 2 Generalized mean
      ('sigma_max             "Maximum singular value")			; ∞ Generalized mean
      ('log_det               "Logarithm of determinant")
      ('entropy               "SVD-based entropy")
      ('storage_efficiency    "Energy storage efficiency")
      ('unobservability_index "Unobservability index")
      ('performance_index     "Performance index"))

    ('system_index "System Index"						;; System index task

      ('cauchy_index        "Cauchy Index")					; Discrete
      ('system_entropy      "System Entropy")					; Linear
      ('gramian_distance    "Gramian Distance")				; Linear
      ('system_symmetry     "System Symmetry")					; Logarithmic
      ('io_coherence        "Input-Output Coherence")				; Logarithmic
      ('system_gain         "System Gain")					; Logarithmic
      ('network_sensitivity "Network Sensitivity")				; Logarithmic
      ('geometric_mean_hsv  "Geometric Mean of Hankel Singular Values")	; Logarithmic
      ('rv_coefficient      "RV Coefficient"))					; Logarithmic

    ('system_norm "Approximate System Norm"					;; System norm task

      ('hardy_inf_norm              "Hardy-Infinity Norm")			; Schatten-1 norm of HSVs
      ('hilbert_schmidt_hankel_norm "Hilbert-Schmidt-Hankel Norm" )		; Schatten-2 norm of HSVs
      ('hankel_norm                 "Hankel Norm")				; Schatten-∞ norm of HSVs
      ('hardy_2_norm                "Hardy-2 Norm"))				; Via Output controllability Gramian
 
    ('tau_function "Tau Function"))						;; Tau function task

  ('Configuration								;;; Configuration Argument

    ('solver "Solver"								;; Integrator for simulation of training and test trajectories

      ('rk1ex   "1st-Order Explicit Runge-Kutta")				; Euler method
      ('rk2ex   "2nd-Order Explicit Runge-Kutta")				; Heun method (DEFAULT)
      ('rk45ex  "4th/5th-Order Explicit Runge-Kutta")				; Dormand-Prince method
      ('@solver "Custom solver"))						; Custom solver function handle with signature @(f,g,t,x0,u,p)

    ('kernel "Inner Product Kernel"						;; Kernel used for the Gramian computation

      ('lie           "Lie Bracket Pseudo Kernel")				; Pseudo kernel for computing Lie bracket of empirical Gramian
      ('hyp           "Hyperbolic-SVD Pseudo Kernel")				; Pseudo kernel for computing hyperbolic SVD difference
      ('sum           "Sum Pseudo Kernel")					; Pseudo kernel for computing the sum of all empirical Gramian elements
      ('trace         "Trace Pseudo Kernel")					; Pseudo kernel for computing the trace of an empirical Gramian
      ('diagonal      "Diagonal Pseudo Kernel")				; Pseudo kernel for computing only the diagonal of an empirical Gramian
      ('position      "Position Pseudo Kernel")				; Pseudo kernel for computing only the upper left diagonal block of half the state-space dimension
      ('velocity      "Velocity Pseudo Kernel")				; Pseudo kernel for computing only the lower right diagonal block of half the state-space dimension
      ('dmd           "DMD Pseudo Kernel")					; Pseudo kernel for computing an approximate Koopman operator via dynamic mode decomposition
      ('linear        "Linear Kernel")						; Linear standard L2 inner product (unit) kernel (DEFAULT)
      ('quadratic     "Quadratic Polynomial Kernel")				; Second order polynomial kernel with unit shift
      ('cubic         "Cubic Polynomial Kernel")				; Third order polynomial kernel with unit shift
      ('sigmoid       "Sigmoid Kernel")					; Sigmoid kernel with unit shift
      ('mercersigmoid "Mercer Sigmoid Kernel")					; Mercer variant of Sigmoid kernel with unit shifts
      ('logarithmic   "Logarithmic Kernel")					; Logarithmic kernel with unit shift
      ('exponential   "Exponential Kernel")					; Exponential kernel
      ('gauss         "Gauss Kernel")						; Gaussian kernel with unit covariance
      ('single        "Single Precision Kernel")				; Linear kernel using single precision floating point numbers
      ('@kernel       "Custom kernel"))					; Custom kernel function handle with signature @(x,y)

    ('training "Training Input"						;; Input used for training trajectories

      ('impulse "Delta Impulse")						; (DEFAULT)
      ('step    "Step Input")							; Constant step input / load vector
      ('chirp   "Decaying Exponential Chirp")					; Exponential chirp with highest frequency 1/dt
      ('sinc    "Cardinal Sine")						; Cardinal sine input
      ('random  "Pseudo-Random Binary"))					; (unseeded) once-sampled scalar {0,1}-sequence

    ('weighting "Trajectory Weighting"						;; Scale or normalize discrete trajectory

      ('none      "No weighting")						; (DEFAULT)
      ('linear    "Linear Time Weighting")					; Scale each trajectory column by the square-root of its time
      ('quadratic "Quadratic Time Weighting")					; Scale each trajectory column by its time
      ('state     "Per-State Weighting")					; Normalize each trajectory column by its 2-norm
      ('scale     "Max-Per-Component Weighting")				; Normalize each state by its inf-norm
      ('rsqrt     "Reciprocal square-root Weighting"))				; Scale each trajectory column by the inverse square-root of its time times pi

    ('centering "Trajectory Centering"						;; Center or shift discrete trajectory

      ('none     "No centering")						; (DEFAULT)
      ('steady   "Steady-State")						; Center (output) trajectory around steady-state (output)
      ('final    "Final State")						; Center (output) trajectory around last state / output
      ('mean     "Arithmetic Mean")						; Center (output) trajectory around temporal arithmetic mean
      ('rms      "Root-Mean-Square")						; Center (output) trajectory around temporal quadratic mean
      ('midrange "Mid-Range"))							; Center (output) trajectory around temporal midrange

    ('scales "Perturbation Scales"						;; Subdivision of perturbation

      ('single      "Single scale")						; (DEFAULT)
      ('linear      "Linear Scaling")						; [0.25,0.50,0.75,1.0] * max_perturbation
      ('geometric   "Geometric Scaling")					; [0.125,0.25,0.5,1.0] * max_perturbation
      ('logarithmic "Logarithmic Scaling")					; [0.001,0.01,0.1,1.0] * max_perturbation
      ('sparse      "Sparse-Grid Scaling"))					; [0.01,0.50,0.99,1.0] * max_perturbation

    ('rotations "Direction Rotations"						;; Perturbation rotations

      ('posneg   "Positive and Negative")					; Unit and negatve unit rotation (DEFAULT)
      ('positive "Only Positive"))						; Only unit rotation

    ('normalization "Gramian Normalization"					;; Improve numerical behavior by scaling

      ('none   "No Normalization")						; (DEFAULT)
      ('steady "Steady-State")							; Normalizes vector field argument and value with steady states
      ('jacobi "Jacobi Preconditioner"))					; Normalizes vector field argument and value with Gramian diagonal

    ('stype "State-Space Gramian Variant"					;; Special Gramian variants

      ('standard "Regular")							; (DEFAULT)
      ('special  "Generic Special")						; Use generic special Gramian variant (for automation purposes)
      ('output_controllability  "Output Controllability Gramian")		; Use with controllability and sensitivity Gramians
      ('average_observability   "Average Observability Gramian")		; Use with observablity and identifiability Gramians
      ('nonsymmetric_minimality "Non-Symmetric Cross Gramian"))		; Use with (linear) cross and joint Gramians

    ('extra_input "Control Explicit Observability"				;; Additional input for observability (initial-state perturbation) simulations

      ('none "No Extra Input")							; # (DEFAULT)
      ('yes  "Use Extra Input"))

    ('pcentering "Parameter Centering"						;; Center parameter and generate parameter perturbations

      ('none        "No Centering")						; Minimal parameter is used as center, scaled linear perturbations (DEFAULT)
      ('linear      "Linear")							; Center at arithmetic mean, linear scaled perturbations
      ('logarithmic "Logarithmic")						; Center at logarithmic mean, logarithmic scaled perturbations (Requires all strictly positive parameter components)
      ('nominal     "Nominal"))						; Center at nominal parameter, linear scaled perturbations

    ('ptype "Parameter-Space Gramian Variant"					;; Special parameter Gramian variants

      ('standard       "Regular")						; Regular parameter Gramian (DEFAULT)
      ('special        "Generic Special")					; Use generic special Gramian variant (for automation purposes)
      ('io_sensitivity "Input-Output Sensitivity Gramian")			; Input-output sensitivity Gramian
      ('coarse_schur   "Coarse Schur Complement")				; Coarse Schur complement computation for empirical identifiability and joint Gramian
      ('exact_schur    "Exact Schur Complement"))				; Exact Schur complement computation for empirical identifiability and joint Gramian

    ('linearity "System Linearity"						;; Assign system linearity

      ('nonlinear "Nonlinear System")						; The system is nonlinear (DEFAULT)
      ('linear    "Linear System")))						; The system is linear, which requires an adjoint vector field
)

