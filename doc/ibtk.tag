<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile>
  <compound kind="class">
    <name>IBTK::AppInitializer</name>
    <filename>class_i_b_t_k_1_1_app_initializer.html</filename>
    <member kind="function">
      <type></type>
      <name>AppInitializer</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>ada9332b85f0d89545aeb328e28080da4</anchor>
      <arglist>(int argc, char *argv[], const std::string &amp;default_log_file_name=&quot;IBAMR.log&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AppInitializer</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>ac3c0ff485bd7f2e7016fe2d84f654fdf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt;</type>
      <name>getInputDatabase</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a93089be6d64ade2c4d670493a61cbb8e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isFromRestart</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a729de7d44c9dda64ea2584a0bb064f4a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt;</type>
      <name>getRestartDatabase</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>ad1b1e95964eb484f85c3e64702e612a5</anchor>
      <arglist>(bool suppress_warning=false)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt;</type>
      <name>getComponentDatabase</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a0e8a40db96caeae108b1f55a46a2a951</anchor>
      <arglist>(const std::string &amp;component_name, bool suppress_warning=false)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dumpVizData</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>ae56d23ab7138d470ab7dab182691fab3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getVizDumpInterval</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>aefb0e8489273cb69e66a4155641ab64d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getVizDumpDirectory</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a74b52f1e3851ce802e2441f09a8ddbe3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::string &gt;</type>
      <name>getVizWriters</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>ac1a12087666f020adf7f2d123c8edaf7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::appu::VisItDataWriter&lt; NDIM &gt; &gt;</type>
      <name>getVisItDataWriter</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>ad53f4de2efe0126f40988d9b6bc25f4b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; LSiloDataWriter &gt;</type>
      <name>getLSiloDataWriter</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>acfc301e6cc246a689473860f5c019b19</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getExodusIIFilename</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a2bd20b743af1a9734f0993f2b409e3a7</anchor>
      <arglist>(const std::string &amp;prefix=&quot;&quot;) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dumpRestartData</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a48df694dea75a8e08392f7e0e5dfd4a5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getRestartDumpInterval</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>ad788770cb84637ed91cd4e0d44fb8068</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getRestartDumpDirectory</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a6e778862026cf9edc61602f007a3b1fb</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dumpPostProcessingData</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a331dba9ac8d3407c94a179b05117545d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getPostProcessingDataDumpInterval</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a8350abf1dcdf5f8d5c79e68cb7b5c9e7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getPostProcessingDataDumpDirectory</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a6d07cb23e2fe05244893e0a0be23c69a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dumpTimerData</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a79a1a371683bf77ec7713ef2afdfaf5c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getTimerDumpInterval</name>
      <anchorfile>class_i_b_t_k_1_1_app_initializer.html</anchorfile>
      <anchor>a5dc228ae54205ab50fba1af30e2e5ca0</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::BGaussSeidelPreconditioner</name>
    <filename>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</filename>
    <base>IBTK::LinearSolver</base>
    <member kind="function">
      <type></type>
      <name>BGaussSeidelPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>a12c9fed6327c847b6ae3710d6ff8b296</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BGaussSeidelPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>aae099fa5167a6cb4be91af2b8411c976</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setComponentPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>a14a9edc888409a97e3da9b56817a2933</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LinearSolver &gt; preconditioner, unsigned int component)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setComponentOperators</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>a79eae016338521e633ee9b8bd020f235</anchor>
      <arglist>(const std::vector&lt; SAMRAI::tbox::Pointer&lt; LinearOperator &gt; &gt; &amp;linear_ops, unsigned int component)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSymmetricPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>a7381d56b2f3c4d17e59965bc2953127c</anchor>
      <arglist>(bool symmetric_preconditioner)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setReversedOrder</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>ad16861526fe07c16792c210018df24fb</anchor>
      <arglist>(bool reverse_order)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>a10712751952deff0d594be9339ae8c3e</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>a622ea4240578328b069f1c67a36c0710</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>a884899c03432ce183ca0eaf611e2a4ab</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setInitialGuessNonzero</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>ac8712e066150e74ad473d8c9f43b9186</anchor>
      <arglist>(bool initial_guess_nonzero=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMaxIterations</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>a764657490107a6d0169d082ad3da8a2a</anchor>
      <arglist>(int max_iterations)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumIterations</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>a0dfe2c117eb697c8de18fc169c94bba4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getResidualNorm</name>
      <anchorfile>class_i_b_t_k_1_1_b_gauss_seidel_preconditioner.html</anchorfile>
      <anchor>a5740421ec01432145e91ba0c29ad5c4a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinearSolver</name>
      <anchorfile>class_i_b_t_k_1_1_linear_solver.html</anchorfile>
      <anchor>a2e22aaea6c75f116c445c7bca446c202</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LinearSolver</name>
      <anchorfile>class_i_b_t_k_1_1_linear_solver.html</anchorfile>
      <anchor>a91367d051f4406590fffa677d6d2a4fe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>getInitialGuessNonzero</name>
      <anchorfile>class_i_b_t_k_1_1_linear_solver.html</anchorfile>
      <anchor>a9106c4895d04ad31893a74525cceb447</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>printClassData</name>
      <anchorfile>class_i_b_t_k_1_1_linear_solver.html</anchorfile>
      <anchor>ab9b7dad18879b274c0efdd9488bc2bda</anchor>
      <arglist>(std::ostream &amp;stream)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setNullspace</name>
      <anchorfile>class_i_b_t_k_1_1_linear_solver.html</anchorfile>
      <anchor>afd6a6761c1e9bc390d0bf7c3177c2538</anchor>
      <arglist>(bool nullspace_contains_constant_vec, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAIVectorReal_NDIM_double &gt; &gt; &amp;nullspace_basis_vecs=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAIVectorReal_NDIM_double &gt; &gt;())</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>getNullspaceContainsConstantVector</name>
      <anchorfile>class_i_b_t_k_1_1_linear_solver.html</anchorfile>
      <anchor>abab55cf1808eba115f5c4784634fc1f4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt; &gt; &amp;</type>
      <name>getNullspaceBasisVectors</name>
      <anchorfile>class_i_b_t_k_1_1_linear_solver.html</anchorfile>
      <anchor>a5cada6e46a790e3009279d2995d71e2c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GeneralSolver</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a22b1dcae6a87c4ca8a56644286a5e3f8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~GeneralSolver</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a12326d0cc29d36764bb79a35b76f4124</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getName</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a0b73449258921a7cafcd3474b42eb386</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>getIsInitialized</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>ae896f1b6a93906574e90f0c5641b0316</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>abf634a801cbe80ffa97519d7f48195dd</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>getHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a231cb072411d956880a6a0b2b511f35d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setSolutionTime</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a833030d1a498cfed5d239ed456a60d25</anchor>
      <arglist>(double solution_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>getSolutionTime</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>af9e70872d61342fc02c1db5c551633a9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>ad99d97a135561487ee01cdc481396a3b</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::pair&lt; double, double &gt;</type>
      <name>getTimeInterval</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a5a6c0ad9bc610ca809254432ecb6ae8c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>getDt</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a76b69886e34948f4945dbc9a04e2cff5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setHierarchyMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a1ebccab5cd5b2f38da57b26bd72e7ed9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; HierarchyMathOps &gt; hier_math_ops)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual SAMRAI::tbox::Pointer&lt; HierarchyMathOps &gt;</type>
      <name>getHierarchyMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a100569d57202bcc4cfdb114c232d855f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getMaxIterations</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a174fc573449b73184d461bef42e29d18</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setAbsoluteTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a36bde17ce3c647eab3fa7952583c3dfd</anchor>
      <arglist>(double abs_residual_tol)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>getAbsoluteTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>aee4b177d922644652955add126881e43</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setRelativeTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a7348a67b3b9c63e3ed6b1f1e0b87ff2a</anchor>
      <arglist>(double rel_residual_tol)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>getRelativeTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a3132938da00894cc5dc83c3dbd7ea83b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setLoggingEnabled</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a95e0355bbc1c512ca1d89be8f80670de</anchor>
      <arglist>(bool enable_logging=true)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>getLoggingEnabled</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a88bcd0fd1f6a2cff1aa98ebcdf00b544</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::BJacobiPreconditioner</name>
    <filename>class_i_b_t_k_1_1_b_jacobi_preconditioner.html</filename>
    <base>IBTK::LinearSolver</base>
    <member kind="function">
      <type></type>
      <name>BJacobiPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_b_jacobi_preconditioner.html</anchorfile>
      <anchor>a181dd9611b9e9d1e4c81a295d77e5bbe</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BJacobiPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_b_jacobi_preconditioner.html</anchorfile>
      <anchor>aaf0827bf4ab69f1a6d3664b269d3402c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setComponentPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_b_jacobi_preconditioner.html</anchorfile>
      <anchor>a33c23a6204be6c7c2967338913aa69b3</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LinearSolver &gt; preconditioner, unsigned int component)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_t_k_1_1_b_jacobi_preconditioner.html</anchorfile>
      <anchor>a544220425f4b4c747cfcd931d34b7977</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_b_jacobi_preconditioner.html</anchorfile>
      <anchor>aafc73e50e5f0097e3d9ee8c4c9337011</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_b_jacobi_preconditioner.html</anchorfile>
      <anchor>a8a4a8efe03b390b00ee78ca35be029fe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setInitialGuessNonzero</name>
      <anchorfile>class_i_b_t_k_1_1_b_jacobi_preconditioner.html</anchorfile>
      <anchor>a939f92d69ed47e9b09b90ed49eb98efd</anchor>
      <arglist>(bool initial_guess_nonzero=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMaxIterations</name>
      <anchorfile>class_i_b_t_k_1_1_b_jacobi_preconditioner.html</anchorfile>
      <anchor>a281ceefa00d3e0ed3621e645b8ccb507</anchor>
      <arglist>(int max_iterations)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumIterations</name>
      <anchorfile>class_i_b_t_k_1_1_b_jacobi_preconditioner.html</anchorfile>
      <anchor>ab69b3ab080d0bfaa32d6082793c1b7d2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getResidualNorm</name>
      <anchorfile>class_i_b_t_k_1_1_b_jacobi_preconditioner.html</anchorfile>
      <anchor>a9344a18cf6d686358290ab583433b853</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartCellDoubleBoundsPreservingConservativeLinearRefine</name>
    <filename>class_i_b_t_k_1_1_cart_cell_double_bounds_preserving_conservative_linear_refine.html</filename>
    <base>RefineOperator&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>CartCellDoubleBoundsPreservingConservativeLinearRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_bounds_preserving_conservative_linear_refine.html</anchorfile>
      <anchor>a92549821a958c4bb7a1768c366126d03</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartCellDoubleBoundsPreservingConservativeLinearRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_bounds_preserving_conservative_linear_refine.html</anchorfile>
      <anchor>a9022919567adee8d3bb96ec9947ae4ef</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>findRefineOperator</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_bounds_preserving_conservative_linear_refine.html</anchorfile>
      <anchor>a1474b4e1773d9fb2b381233507c2abdd</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;var, const std::string &amp;op_name) const </arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getOperatorName</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_bounds_preserving_conservative_linear_refine.html</anchorfile>
      <anchor>a9ee469780a78bf8d9604c9fc0056c024</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getOperatorPriority</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_bounds_preserving_conservative_linear_refine.html</anchorfile>
      <anchor>a69ee64471ee022a84e64eba5e03f2bdc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_bounds_preserving_conservative_linear_refine.html</anchorfile>
      <anchor>a432f65e9ace5801b6f7754388a7f973d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>refine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_bounds_preserving_conservative_linear_refine.html</anchorfile>
      <anchor>ab55716da98cbfa967e77507e54b105b9</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, int dst_component, int src_component, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartCellDoubleCubicCoarsen</name>
    <filename>class_i_b_t_k_1_1_cart_cell_double_cubic_coarsen.html</filename>
    <base>CoarsenOperator&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>CartCellDoubleCubicCoarsen</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_cubic_coarsen.html</anchorfile>
      <anchor>a0ca9f2ce4ed49b6ff8a0cc3250501f58</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartCellDoubleCubicCoarsen</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_cubic_coarsen.html</anchorfile>
      <anchor>ad6ea23cd10074358a96a7067c740a220</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>findCoarsenOperator</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_cubic_coarsen.html</anchorfile>
      <anchor>a95333937aa6bd67f064a73b68722c6d6</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;var, const std::string &amp;op_name) const </arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getOperatorName</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_cubic_coarsen.html</anchorfile>
      <anchor>a526f97e16b5a87bda7ff360da2c119c2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getOperatorPriority</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_cubic_coarsen.html</anchorfile>
      <anchor>adcf0c8a332de71f3896320b741c90599</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_cubic_coarsen.html</anchorfile>
      <anchor>ab4410b1992736a6f7a5f93a529dd1740</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>coarsen</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_cubic_coarsen.html</anchorfile>
      <anchor>ac24c78176dae0b69eff5e5cb2c56a657</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, int dst_component, int src_component, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;coarse_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartCellDoubleLinearCFInterpolation</name>
    <filename>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</filename>
    <base>IBTK::CoarseFineBoundaryRefinePatchStrategy</base>
    <member kind="function">
      <type></type>
      <name>CartCellDoubleLinearCFInterpolation</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>a64afb80c89c78bac3ded5be9ab86e18f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartCellDoubleLinearCFInterpolation</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>a845b04dd078f0525cee78eb3f6c291cb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>a6f7c2ad124f720b8fe908f8110a8d571</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getRefineOpStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>a17f5535b1dc6d674ea690babf9ff0be9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>aef874e7f6181a15858c828bfc1ba1331</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>aad3736f35a5afc7d8ec730836f5a4f81</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConsistentInterpolationScheme</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>ab687379f1aaaed52b18be59a6f584a51</anchor>
      <arglist>(bool consistent_type_2_bdry)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndex</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>a04068df7c8233eae5c571809aaccade7</anchor>
      <arglist>(int patch_data_index)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>a3ee4c7f70a8eab933df0d17eed306ef1</anchor>
      <arglist>(const std::set&lt; int &gt; &amp;patch_data_indices)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>a0bf56c7f90117bdc6ce6d41de2b45794</anchor>
      <arglist>(const SAMRAI::hier::ComponentSelector &amp;patch_data_indices)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>ac8e1d32656b9faa236856506e56ba5fa</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>a4d034ceb130b383e999feb458fc62121</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeNormalExtension</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_linear_c_f_interpolation.html</anchorfile>
      <anchor>aa5251c8ed3ae775d52746d847023d983</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CoarseFineBoundaryRefinePatchStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>ab26ae1d7544a3f051e2c68c8d58a0fc5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~CoarseFineBoundaryRefinePatchStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>a8ce9953beb04ff8089e881415d49a14e</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartCellDoubleQuadraticCFInterpolation</name>
    <filename>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</filename>
    <base>IBTK::CoarseFineBoundaryRefinePatchStrategy</base>
    <member kind="function">
      <type></type>
      <name>CartCellDoubleQuadraticCFInterpolation</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a681e32ded9dd6e3782657da21c27ae60</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartCellDoubleQuadraticCFInterpolation</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a0d0628c76e04d0a2b4c2e5f754a8660f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>aedfcc63441e4c4e3c2316300ef85f57b</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getRefineOpStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a7eb889b7e2a45f3e594bd0ddc3b85a41</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a33b40f02751a109a6c10a0f4da7e2814</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a592949298148f48a34eebf376d7b946e</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConsistentInterpolationScheme</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a38cabcf15dc8b878a08e2ac1403550b5</anchor>
      <arglist>(bool consistent_type_2_bdry)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndex</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a5384c98ee7279351da27c1ec323a083c</anchor>
      <arglist>(int patch_data_index)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a0aa56a98a4be6164546e10a76c62e011</anchor>
      <arglist>(const std::set&lt; int &gt; &amp;patch_data_indices)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>ade9279a2344ad593d89e959fa8b95205</anchor>
      <arglist>(const SAMRAI::hier::ComponentSelector &amp;patch_data_indices)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a7a7195c2f7b04fe6a7feda866a02e6cb</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a102f0660175d0372f949704de55aa5b4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeNormalExtension</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a35fa7f45851cd30096807f23f5cc0cfb</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartCellDoubleQuadraticRefine</name>
    <filename>class_i_b_t_k_1_1_cart_cell_double_quadratic_refine.html</filename>
    <base>RefineOperator&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>CartCellDoubleQuadraticRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_refine.html</anchorfile>
      <anchor>a2d049a7640483b1f25d2f1c243a6d304</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartCellDoubleQuadraticRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_refine.html</anchorfile>
      <anchor>a3a636bf4a9919bcdb06065291999f1ec</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>findRefineOperator</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_refine.html</anchorfile>
      <anchor>aef5bdcd3982cd570a6de20fcdc37ed39</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;var, const std::string &amp;op_name) const </arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getOperatorName</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_refine.html</anchorfile>
      <anchor>af5496e34910949aba641770d5110fc51</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getOperatorPriority</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_refine.html</anchorfile>
      <anchor>a00a6f656be9b3b7ba4d94afb66417f53</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_refine.html</anchorfile>
      <anchor>ace6c3070c53243275c5b0a7d33387b52</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>refine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_double_quadratic_refine.html</anchorfile>
      <anchor>a521a763ce7966b78d4651b4ad7a384f6</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, int dst_component, int src_component, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartCellRobinPhysBdryOp</name>
    <filename>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</filename>
    <base>IBTK::RobinPhysBdryPatchStrategy</base>
    <member kind="function">
      <type></type>
      <name>CartCellRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</anchorfile>
      <anchor>a25a062ffeb191e47cc65e91f896ff7c8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartCellRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</anchorfile>
      <anchor>a7a9e8e240923a054749b6f326420d02c</anchor>
      <arglist>(int patch_data_index, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartCellRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</anchorfile>
      <anchor>aee505bb57b7033a006a57f0abe31af35</anchor>
      <arglist>(const std::set&lt; int &gt; &amp;patch_data_indices, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartCellRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</anchorfile>
      <anchor>a60dd12549d1fe5bf60ade9425d214a25</anchor>
      <arglist>(const SAMRAI::hier::ComponentSelector &amp;patch_data_indices, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartCellRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</anchorfile>
      <anchor>aaa6eae31bb638c4f506c27c5b524e874</anchor>
      <arglist>(int patch_data_index, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartCellRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</anchorfile>
      <anchor>a9a5c205830d654c393ad337f608a3fa8</anchor>
      <arglist>(const std::set&lt; int &gt; &amp;patch_data_indices, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartCellRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</anchorfile>
      <anchor>aa0e618941fd7649237f4677f707772d1</anchor>
      <arglist>(const SAMRAI::hier::ComponentSelector &amp;patch_data_indices, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartCellRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</anchorfile>
      <anchor>a5d636f5f1452ec54e063661db58affbd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>accumulateFromPhysicalBoundaryData</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</anchorfile>
      <anchor>a1fd4cf1faee55731d4cc9623a22525dc</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</anchorfile>
      <anchor>a68ef4841bc16ad4ea578753b09426284</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getRefineOpStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_cell_robin_phys_bdry_op.html</anchorfile>
      <anchor>acce72fd5d6ef05aa0e07835527ffada6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RobinPhysBdryPatchStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>a1eca916c8e5aee9bc90271ee139de29f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RobinPhysBdryPatchStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>a84b7296feb9026cb3716a7a1bd58cc1b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndex</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>a2d0e8fa290f85719e7f5099c01b63773</anchor>
      <arglist>(int patch_data_index)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>a052f712ac937e1c7a571530a194b68e4</anchor>
      <arglist>(const std::set&lt; int &gt; &amp;patch_data_indices)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>a20127fe01e0dfac65464429a5dc42614</anchor>
      <arglist>(const SAMRAI::hier::ComponentSelector &amp;patch_data_indices)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoef</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>ae34cba46245c65d607c8949ae683d8cb</anchor>
      <arglist>(SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>aac101c68fa7eabdc08f70016895fae22</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>a0ccb94ff10c87a6b2099732f6b1c8761</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>getHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>a365f10451425768e5c316633a46cfb97</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>ae9e98f8407efd77fd69d8374e0a672ed</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>a88dd5fd0097e7698b23d71f9335a1030</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartExtrapPhysBdryOp</name>
    <filename>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</filename>
    <base>RefinePatchStrategy&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>CartExtrapPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>a73447094b74acf6eeeb02e4817a533f1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartExtrapPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>a1e9e62d5e61cb51929edd375a7931f06</anchor>
      <arglist>(int patch_data_index, const std::string &amp;extrap_type=&quot;CONSTANT&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartExtrapPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>a88a73209d269cab1c3b3fae781c4abb8</anchor>
      <arglist>(const std::set&lt; int &gt; &amp;patch_data_indices, const std::string &amp;extrap_type=&quot;CONSTANT&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartExtrapPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>a9d76ee5ca83cfe441f5e134976ac913d</anchor>
      <arglist>(const SAMRAI::hier::ComponentSelector &amp;patch_data_indices, const std::string &amp;extrap_type=&quot;CONSTANT&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartExtrapPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>a9c9e6036b4760c5b52f46509cc2cc89e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndex</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>a5e893d95f5bc4e8107d6c6ac9576408f</anchor>
      <arglist>(int patch_data_index)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>a6a876257eacec244f18c4c35b8598c58</anchor>
      <arglist>(const std::set&lt; int &gt; &amp;patch_data_indices)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>a4423957c229bd09135cad6d3fafb509e</anchor>
      <arglist>(const SAMRAI::hier::ComponentSelector &amp;patch_data_indices)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setExtrapolationType</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>ac7b093fad21f21288ac7d894329b54fc</anchor>
      <arglist>(const std::string &amp;extrap_type)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>a23dfb1943f96cf2993e321cfceaf51d1</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getRefineOpStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>a0a58d7d341bb6830afe92195c50d6826</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>ac864a80889867f9bb7930b60bba083e2</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_extrap_phys_bdry_op.html</anchorfile>
      <anchor>afb8abff751631ea77321675969aee171</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartGridFunction</name>
    <filename>class_i_b_t_k_1_1_cart_grid_function.html</filename>
    <member kind="function">
      <type></type>
      <name>CartGridFunction</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function.html</anchorfile>
      <anchor>ab0f454d12fbaa3d0b1ed3484407dc288</anchor>
      <arglist>(const std::string &amp;object_name=&quot;&quot;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~CartGridFunction</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function.html</anchorfile>
      <anchor>ab54a35b6123ddbc00fe124eac49efb22</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>isTimeDependent</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function.html</anchorfile>
      <anchor>aa3aebc3754f2b461964e4b38b00fc68c</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setDataOnPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function.html</anchorfile>
      <anchor>a0e668b2f6dab2bde34bb02a76cc4908c</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, double data_time, bool initial_time=false, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setDataOnPatchLevel</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function.html</anchorfile>
      <anchor>a7cfb06ca1811242db0bea45d75ceea43</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level, double data_time, bool initial_time=false)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setDataOnPatch</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function.html</anchorfile>
      <anchor>afc47de5fa45d121ce21550a25fcd6cb4</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, double data_time, bool initial_time=false, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt;(NULL))=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartGridFunctionSet</name>
    <filename>class_i_b_t_k_1_1_cart_grid_function_set.html</filename>
    <base>IBTK::CartGridFunction</base>
    <member kind="function">
      <type></type>
      <name>CartGridFunctionSet</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function_set.html</anchorfile>
      <anchor>aca4fa3b7c4c39f4dc1476ec75a097a60</anchor>
      <arglist>(const std::string &amp;object_name=&quot;&quot;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~CartGridFunctionSet</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function_set.html</anchorfile>
      <anchor>a9a6695bcb8a58832339866fe2c83d17a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>addFunction</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function_set.html</anchorfile>
      <anchor>af6f6cd4cc378202e46bc0776ce7380d1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; CartGridFunction &gt; fcn)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isTimeDependent</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function_set.html</anchorfile>
      <anchor>a4523981cdc45f7896b33462be59b0487</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function_set.html</anchorfile>
      <anchor>aa4c561029d7227ed16feb7a45b4ef581</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, double data_time, bool initial_time=false, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatchLevel</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function_set.html</anchorfile>
      <anchor>adb36db29fbd2e290fca0ba661c9c5c40</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level, double data_time, bool initial_time=false)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatch</name>
      <anchorfile>class_i_b_t_k_1_1_cart_grid_function_set.html</anchorfile>
      <anchor>aaa8622e9a02186c89089e1d310eaaff7</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, double data_time, bool initial_time=false, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt;(NULL))</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartSideDoubleCubicCoarsen</name>
    <filename>class_i_b_t_k_1_1_cart_side_double_cubic_coarsen.html</filename>
    <base>CoarsenOperator&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>CartSideDoubleCubicCoarsen</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_cubic_coarsen.html</anchorfile>
      <anchor>a3e8e69321eaae62a3191588b83ee77b6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartSideDoubleCubicCoarsen</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_cubic_coarsen.html</anchorfile>
      <anchor>a1188d42c684c3e656573af3b9cbebcb9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>findCoarsenOperator</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_cubic_coarsen.html</anchorfile>
      <anchor>a733f97618955c171d3408b0b9797ed52</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;var, const std::string &amp;op_name) const </arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getOperatorName</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_cubic_coarsen.html</anchorfile>
      <anchor>a96e19fb5261e051b3f84cd22670bc977</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getOperatorPriority</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_cubic_coarsen.html</anchorfile>
      <anchor>ab3b05d43ea002d7cbdc3b57ca9a706df</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_cubic_coarsen.html</anchorfile>
      <anchor>a8bdce96163d553f5dec7856a5ba03e5c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>coarsen</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_cubic_coarsen.html</anchorfile>
      <anchor>ad39230ef851880a7a15e2129ebe6eb7b</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, int dst_component, int src_component, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;coarse_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartSideDoubleDivPreservingRefine</name>
    <filename>class_i_b_t_k_1_1_cart_side_double_div_preserving_refine.html</filename>
    <base>RefinePatchStrategy&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>CartSideDoubleDivPreservingRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_div_preserving_refine.html</anchorfile>
      <anchor>aa38ee4c97182baad5931064de3489963</anchor>
      <arglist>(int u_dst_idx, int u_src_idx, int indicator_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineOperator&lt; NDIM &gt; &gt; refine_op, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenOperator&lt; NDIM &gt; &gt; coarsen_op, double fill_time, SAMRAI::xfer::RefinePatchStrategy&lt; NDIM &gt; *phys_bdry_op)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~CartSideDoubleDivPreservingRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_div_preserving_refine.html</anchorfile>
      <anchor>afbe3520a3978ef5320083f4be12c4a04</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_div_preserving_refine.html</anchorfile>
      <anchor>a713e5d16ddeedf8ab3d05486172406f0</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getRefineOpStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_div_preserving_refine.html</anchorfile>
      <anchor>a4d87bd55aaca946792f025e9a224447d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>preprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_div_preserving_refine.html</anchorfile>
      <anchor>a361fc0604e9ef82ab8487b2c573e89fe</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>postprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_div_preserving_refine.html</anchorfile>
      <anchor>ac7c3086eecf2929ab7345f49873f2059</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>REFINE_OP_STENCIL_WIDTH</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_div_preserving_refine.html</anchorfile>
      <anchor>abbbdef1d951213bb759e9884fcb89113</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartSideDoubleQuadraticCFInterpolation</name>
    <filename>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</filename>
    <base>IBTK::CoarseFineBoundaryRefinePatchStrategy</base>
    <member kind="function">
      <type></type>
      <name>CartSideDoubleQuadraticCFInterpolation</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a34f7b635e6204a6ef95768df62933067</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartSideDoubleQuadraticCFInterpolation</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a49b9b141733d641428b94704d0b28c5a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a12f420ae10635d0da8e09fd1f7462fcd</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getRefineOpStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>af730b584430b0e47392aebec5f9f8194</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>aa33220f8d4f13ff8e255e701eaae37de</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a6fc4864b72add7d3d0f4ecc43af5b394</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setConsistentInterpolationScheme</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a0a5a37e56f510bc2adb4d641836d62ff</anchor>
      <arglist>(bool consistent_type_2_bdry)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndex</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a8e19aa13e6f7e29ba37489188e32bcb9</anchor>
      <arglist>(int patch_data_index)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>ab6ad41c3eaecb7f25edc8cb23fbd180d</anchor>
      <arglist>(const std::set&lt; int &gt; &amp;patch_data_indices)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a23c616862bf3b921583e58d54beb9edb</anchor>
      <arglist>(const SAMRAI::hier::ComponentSelector &amp;patch_data_indices)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>aafa6d19f3537f1a99ffbbcd2632d78dd</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clearPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a6fd4d6e2f287f2cf33110242666a51bf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeNormalExtension</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_quadratic_c_f_interpolation.html</anchorfile>
      <anchor>a9451a13ff248d4e10ec425a3da27ae71</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartSideDoubleSpecializedConstantRefine</name>
    <filename>class_i_b_t_k_1_1_cart_side_double_specialized_constant_refine.html</filename>
    <base>RefineOperator&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>CartSideDoubleSpecializedConstantRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_constant_refine.html</anchorfile>
      <anchor>afe5a9c700c25fd08b02c789859ea771f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartSideDoubleSpecializedConstantRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_constant_refine.html</anchorfile>
      <anchor>ad04a7ddaab4384210f52c65c493d4468</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>findRefineOperator</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_constant_refine.html</anchorfile>
      <anchor>a0bb2c2228663a2a7cfb99eaa5f94f2bb</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;var, const std::string &amp;op_name) const </arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getOperatorName</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_constant_refine.html</anchorfile>
      <anchor>a0933a859431ae7de02c3a560e5845e76</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getOperatorPriority</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_constant_refine.html</anchorfile>
      <anchor>a440aa6c0097aeb41c99cd6f85affd66e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_constant_refine.html</anchorfile>
      <anchor>a8e6af190edd5b15f0c120352808d6e2d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>refine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_constant_refine.html</anchorfile>
      <anchor>ace4fafb0989ca618c83f3dab127837dd</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, int dst_component, int src_component, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartSideDoubleSpecializedLinearRefine</name>
    <filename>class_i_b_t_k_1_1_cart_side_double_specialized_linear_refine.html</filename>
    <base>RefineOperator&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>CartSideDoubleSpecializedLinearRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_linear_refine.html</anchorfile>
      <anchor>a7921196670f48b3c254c1cf92a3b84f8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartSideDoubleSpecializedLinearRefine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_linear_refine.html</anchorfile>
      <anchor>abf773bb5bb244a5967ce36f7c713a024</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>findRefineOperator</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_linear_refine.html</anchorfile>
      <anchor>a1f74eab14f9ab97b0b4ddeb864707888</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;var, const std::string &amp;op_name) const </arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getOperatorName</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_linear_refine.html</anchorfile>
      <anchor>acc82df306592ba7ca2a528ada19ba924</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getOperatorPriority</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_linear_refine.html</anchorfile>
      <anchor>aac12da07397fb6dca6a15b88fb3f8cba</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_linear_refine.html</anchorfile>
      <anchor>aad931662a0bb8b0ab4355cf29a27f8e6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>refine</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_double_specialized_linear_refine.html</anchorfile>
      <anchor>afdc862f5364f1aaa8ab5ef3b0b933d17</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, int dst_component, int src_component, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CartSideRobinPhysBdryOp</name>
    <filename>class_i_b_t_k_1_1_cart_side_robin_phys_bdry_op.html</filename>
    <base>IBTK::RobinPhysBdryPatchStrategy</base>
    <member kind="function">
      <type></type>
      <name>CartSideRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_robin_phys_bdry_op.html</anchorfile>
      <anchor>ae9a5978cd4b785cc2b36edc10a4ecaa0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartSideRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_robin_phys_bdry_op.html</anchorfile>
      <anchor>adcbfa1989e13a10e7d579e20e89f5766</anchor>
      <arglist>(int patch_data_index, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartSideRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_robin_phys_bdry_op.html</anchorfile>
      <anchor>ac9f25a7c0ba3a8dc16ed7224ff7a1771</anchor>
      <arglist>(const std::set&lt; int &gt; &amp;patch_data_indices, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartSideRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_robin_phys_bdry_op.html</anchorfile>
      <anchor>a1f2be68403978167eb4b612f6ec9341e</anchor>
      <arglist>(const SAMRAI::hier::ComponentSelector &amp;patch_data_indices, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CartSideRobinPhysBdryOp</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_robin_phys_bdry_op.html</anchorfile>
      <anchor>aa0397f0285e8537f88bf61d8f3480a11</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>accumulateFromPhysicalBoundaryData</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_robin_phys_bdry_op.html</anchorfile>
      <anchor>acfe94e50242d16c838e2db9bdec761c1</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_robin_phys_bdry_op.html</anchorfile>
      <anchor>a2816be47b9897bb02d81e59112ad032c</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getRefineOpStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cart_side_robin_phys_bdry_op.html</anchorfile>
      <anchor>a508ba0828644d2c702c98ce4fa7723f6</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CCLaplaceOperator</name>
    <filename>class_i_b_t_k_1_1_c_c_laplace_operator.html</filename>
    <base>IBTK::LaplaceOperator</base>
    <member kind="function">
      <type></type>
      <name>CCLaplaceOperator</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_laplace_operator.html</anchorfile>
      <anchor>a7447536881aafee6eca9b9b70042d1ed</anchor>
      <arglist>(const std::string &amp;object_name, bool homogeneous_bc=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CCLaplaceOperator</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_laplace_operator.html</anchorfile>
      <anchor>a91c5556bb785b7145299c79da3fc2c59</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>apply</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_laplace_operator.html</anchorfile>
      <anchor>ab4c1bba4174b9354fd000d2fc26005bb</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_laplace_operator.html</anchorfile>
      <anchor>a44eb144df1e7b32b645ffcb9b86c1183</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_laplace_operator.html</anchorfile>
      <anchor>ab2bf1837e364b5d75ff98361df313e58</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LaplaceOperator</name>
      <anchorfile>class_i_b_t_k_1_1_laplace_operator.html</anchorfile>
      <anchor>aeda12d8de2c7bd1a3037126d759438af</anchor>
      <arglist>(const std::string &amp;object_name, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LaplaceOperator</name>
      <anchorfile>class_i_b_t_k_1_1_laplace_operator.html</anchorfile>
      <anchor>a5c8f5824d3be7cbbdb37dd18379b8bd2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPoissonSpecifications</name>
      <anchorfile>class_i_b_t_k_1_1_laplace_operator.html</anchorfile>
      <anchor>a0c19d8cd5a8c3fca266994ca66bd14e6</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const SAMRAI::solv::PoissonSpecifications &amp;</type>
      <name>getPoissonSpecifications</name>
      <anchorfile>class_i_b_t_k_1_1_laplace_operator.html</anchorfile>
      <anchor>a4c25ce7f38dbe97a75c42c3830bd10dc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBcCoef</name>
      <anchorfile>class_i_b_t_k_1_1_laplace_operator.html</anchorfile>
      <anchor>a436d9612f80a69e5c0e16e4ef96d7f78</anchor>
      <arglist>(SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_t_k_1_1_laplace_operator.html</anchorfile>
      <anchor>a39c608f4381f532d763df9508a341449</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;</type>
      <name>getPhysicalBcCoefs</name>
      <anchorfile>class_i_b_t_k_1_1_laplace_operator.html</anchorfile>
      <anchor>ab05912d4717ea62e5b68e58e5f14e84c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinearOperator</name>
      <anchorfile>class_i_b_t_k_1_1_linear_operator.html</anchorfile>
      <anchor>ad06817b99a95cf3e3719c6156f4c472f</anchor>
      <arglist>(const std::string &amp;object_name, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LinearOperator</name>
      <anchorfile>class_i_b_t_k_1_1_linear_operator.html</anchorfile>
      <anchor>ad7f1730a4aa2346452a1707ddeea393e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>modifyRhsForInhomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_linear_operator.html</anchorfile>
      <anchor>a3be87569ebc5fd75763737a7cf96d8f7</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GeneralOperator</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a7de0bee126a36813ff46acdd35a83fbd</anchor>
      <arglist>(const std::string &amp;object_name, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~GeneralOperator</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>adabaa62c9e26e23e8cd1b913edaec917</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getName</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a359574a653b1d298382ac026118fd4fa</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>getIsInitialized</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a75ac38292d2822b6009e2f26fe8a4213</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>aeacf5e7b22d41aedf019bcbf120d4eec</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>getHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a8cf1da9e3c3ec5b8b4f5b4bd138ecb3f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setSolutionTime</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a88abf2c68b224d8a0f91d6d299b75d21</anchor>
      <arglist>(double solution_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>getSolutionTime</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a806432dfd8c0cb278d441fc8d9157c4c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a5f516bfb19305218d943ceb0baca8a84</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::pair&lt; double, double &gt;</type>
      <name>getTimeInterval</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a19e4571413593e2fc06adb694b1793cf</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>getDt</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>afa759cf327504e0bdad8f4bd3511c33f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setHierarchyMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a780ac0d666e30e62d2700deed54751ac</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; HierarchyMathOps &gt; hier_math_ops)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual SAMRAI::tbox::Pointer&lt; HierarchyMathOps &gt;</type>
      <name>getHierarchyMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>ad7377a5ebfe0765d50facf6b10f0464f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>applyAdd</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a47ab21cc0e4e223257497477bb9dec66</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;z)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setLoggingEnabled</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a2335f8ca447b1dd55d763aca7d6d2ee4</anchor>
      <arglist>(bool enable_logging=true)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>getLoggingEnabled</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a74e9130f71f0364d8f4b00b30583027e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>printClassData</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>ae3eab864386dd598c5d33563d21078d9</anchor>
      <arglist>(std::ostream &amp;stream)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CCPoissonHypreLevelSolver</name>
    <filename>class_i_b_t_k_1_1_c_c_poisson_hypre_level_solver.html</filename>
    <base>IBTK::LinearSolver</base>
    <base>IBTK::PoissonSolver</base>
    <member kind="function">
      <type></type>
      <name>CCPoissonHypreLevelSolver</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>aa9a07e8fd6912e6df9fb1c9f22c5145f</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CCPoissonHypreLevelSolver</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>a85ebf5bcf82d390fa5cff6d9239e6ca4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>a48468b20fd51c5284853fb93c83761dd</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>a9e511ad3f4dc2919be7e58d5a7659a5d</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>a553a6960ca92ed5f543eee6ddb352148</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setInitialGuessNonzero</name>
      <anchorfile>class_i_b_t_k_1_1_linear_solver.html</anchorfile>
      <anchor>aaaa23c550ec4ef095c9847e6998086bf</anchor>
      <arglist>(bool initial_guess_nonzero=true)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setMaxIterations</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a8b1b84a979601c15789ea4ed32a280d8</anchor>
      <arglist>(int max_iterations)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getNumIterations</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>aefd17fcb2d0f1e6f455a3fcfea0f9acc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>getResidualNorm</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a6320020fcde8abeeb795443bc07a5940</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PoissonSolver</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_solver.html</anchorfile>
      <anchor>a0ea307927dc2b710f037d7f9d07ae055</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PoissonSolver</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_solver.html</anchorfile>
      <anchor>ab24d1148b73f05634b2a9e6f808f9548</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPoissonSpecifications</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_solver.html</anchorfile>
      <anchor>aea712d9312d7ba336ff0e38b4647dd73</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBcCoef</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_solver.html</anchorfile>
      <anchor>a1723a1776500fcb813c0b956f74ed9b2</anchor>
      <arglist>(SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_solver.html</anchorfile>
      <anchor>a5fb0da25906e0102bc919db74e07b40d</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;</type>
      <name>allocate_solver</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>a090276e3c7a7ece8fc6625d6068d8227</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CCPoissonPETScLevelSolver</name>
    <filename>class_i_b_t_k_1_1_c_c_poisson_p_e_t_sc_level_solver.html</filename>
    <base>IBTK::PETScLevelSolver</base>
    <base>IBTK::PoissonSolver</base>
    <member kind="function">
      <type></type>
      <name>CCPoissonPETScLevelSolver</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a9b48c1a0e95350c2eac97719d58fbf87</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CCPoissonPETScLevelSolver</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a82fd858233842cbae27b6d15d0eadc96</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PETScLevelSolver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>afa84819b952a8f270e8995976437225b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PETScLevelSolver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a9aedb724502c0d5362fccf214bcf5f6f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setKSPType</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a2f17c2f0d0ccadd6dce099de97a9ec50</anchor>
      <arglist>(const std::string &amp;ksp_type)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setOptionsPrefix</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a8afe46c4cae67979b8f2ec6837959424</anchor>
      <arglist>(const std::string &amp;options_prefix)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setNullspace</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a789ffba3e577eaa0db283a01dbf7290c</anchor>
      <arglist>(bool contains_constant_vec, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAIVectorReal_NDIM_double &gt; &gt; &amp;nullspace_basis_vecs=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAIVectorReal_NDIM_double &gt; &gt;())</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a4e28e8215c8693a0c0f7a6483d22c9d0</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>af336e4fe08f6b0182520042357d96838</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a0d3e5807a69529dfb2bf208f35c1f4bf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;</type>
      <name>allocate_solver</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a1126410dd077534831b72cae7d9d56c6</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initializeSolverStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>ab90ac94ee3f0af29d0cfe62109f9ccc4</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>deallocateSolverStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a4563811e8359c315e142cf6d2b147596</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>copyToPETScVec</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a756693440eae92b106569c3b3b6322df</anchor>
      <arglist>(Vec &amp;petsc_x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>copyFromPETScVec</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a0b35b1299d009d534109292dd0d49a85</anchor>
      <arglist>(Vec &amp;petsc_x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>setupKSPVecs</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a81fa2b52379db705376d2ef5a0a7fac5</anchor>
      <arglist>(Vec &amp;petsc_x, Vec &amp;petsc_b, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>init</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>aad661c41bf6480c2bf712c40f2d236fc</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>setupNullspace</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a867a0069eec77ef33407d3ff3abfc4b4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt;</type>
      <name>d_hierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a3198c09483925129caa99c8cb5aed0a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>d_level_num</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a928a175bddd45f848df594021ba5fbe5</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CCPoissonPointRelaxationFACOperator</name>
    <filename>class_i_b_t_k_1_1_c_c_poisson_point_relaxation_f_a_c_operator.html</filename>
    <base>IBTK::PoissonFACPreconditionerStrategy</base>
    <member kind="function">
      <type></type>
      <name>CCPoissonPointRelaxationFACOperator</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>adec655b90da51d9c694bc245e4f5c17b</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CCPoissonPointRelaxationFACOperator</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a8d165f6817e748b16341f9bd39bf708f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSmootherType</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a497cb757f32d4c763a5f788b8a0864aa</anchor>
      <arglist>(const std::string &amp;smoother_type)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCoarseSolverType</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a725a5a96c3d5c72c5f8c3f19e15f4c0b</anchor>
      <arglist>(const std::string &amp;coarse_solver_type)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>smoothError</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a0dd50131a02d8bb36ac134188d9aa1de</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;error, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, int level_num, int num_sweeps, bool performing_pre_sweeps, bool performing_post_sweeps)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveCoarsestLevel</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>acbc5078cbb3dce51f87be0a589dca91a</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;error, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, int coarsest_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeResidual</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a32dbd666694bc5f726a1a6408241e083</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs, int coarsest_level_num, int finest_level_num)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PoissonFACPreconditionerStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a8b75d4dde70aa3e92ce04f891e281322</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; scratch_var, int ghost_cell_width, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PoissonFACPreconditionerStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a643beb19c9d9dc2ece18e6ff4321dbc2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPoissonSpecifications</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>af19a4096ad3764c9b1861ef5654ed9f7</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBcCoef</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ad66361b76a2dbc8bfc02198df765417a</anchor>
      <arglist>(SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a76d64d123750044c2f1ee35ead41b148</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setResetLevels</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a4f50e2f12015aa107ff470d0cc759af9</anchor>
      <arglist>(int coarsest_ln, int finest_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCoarseSolverMaxIterations</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ab0f5b2ac0d2dbfbd7714dacdbe941bbb</anchor>
      <arglist>(int coarse_solver_max_iterations)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCoarseSolverAbsoluteTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a99149b364639fd308002549627851479</anchor>
      <arglist>(double coarse_solver_abs_residual_tol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCoarseSolverRelativeTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>af43929a3199f528d4f9cfc9462d474c4</anchor>
      <arglist>(double coarse_solver_rel_residual_tol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setProlongationMethod</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a867908df9d9ed3c7b7dd9fdf548982da</anchor>
      <arglist>(const std::string &amp;prolongation_method)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setRestrictionMethod</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a374ac35eb38f9568d279aef4dcec21bd</anchor>
      <arglist>(const std::string &amp;restriction_method)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>restrictResidual</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>aec8b4cfe6e2fcaa81711c82c676301ca</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;src, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;dst, int dst_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>prolongError</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a855133baaadfc65697765758da75e4aa</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;src, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;dst, int dst_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>prolongErrorAndCorrect</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a6b4620aa5a96be68073b030f8cc743a1</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;src, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;dst, int dst_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a1f8e6dc9e684e4e49ab91313ef169df4</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>af457fa0e3c94a97609db59f9c88e6c9f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>allocateScratchData</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a2d4b4f3e9b9d5d8c76e7afb3dd3c9252</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateScratchData</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>aed54068975f32aa66f20f1289bce3718</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FACPreconditionerStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a36b00bb2e0f29b521ce8625bcc0ad56f</anchor>
      <arglist>(const std::string &amp;object_name, bool homogeneous_bc=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~FACPreconditionerStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ac072ab69fd4e87f6ffbf7fe08180e625</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getName</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a87be932de79c3fce33f0e6b5bfe8a3f8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>getIsInitialized</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>acbb26d6f7158fc73b115f52fe2815a23</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setFACPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a8594d349067f96c4054c74f59202d642</anchor>
      <arglist>(SAMRAI::tbox::ConstPointer&lt; FACPreconditioner &gt; preconditioner)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a5179b59b8dd954bb94a5646539c5277f</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>getHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ac797a92ff4d63df7df9652ea9f7b2754</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setSolutionTime</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a6d780177aac907edac2486ed23a25ac4</anchor>
      <arglist>(double solution_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>getSolutionTime</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a42d579dc6f016cf3b380d804f6146145</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a5884b78f7d62010737d8bcdd55f4f6de</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::pair&lt; double, double &gt;</type>
      <name>getTimeInterval</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a286d5a403980cb599263d4d94a2685fe</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>getDt</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a145619216f4d0456e050d8a9eb7adda6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>printClassData</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>aef8f224c8b123e0c9c7637c404c694ba</anchor>
      <arglist>(std::ostream &amp;stream)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;</type>
      <name>allocate_solver</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a6e4d320f84001d6f8880638de3034302</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initializeOperatorStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a07cada8bf42294a38651b97cb022d482</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs, int coarsest_reset_ln, int finest_reset_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>deallocateOperatorStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>acce2d5e9c490425e7dae19f56401029a</anchor>
      <arglist>(int coarsest_reset_ln, int finest_reset_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>xeqScheduleProlongation</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a11310a5306d778f58d69cfcfc0d7338e</anchor>
      <arglist>(int dst_idx, int src_idx, int dst_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>xeqScheduleRestriction</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a02869d1a2ea3e5b75091b84344683c00</anchor>
      <arglist>(int dst_idx, int src_idx, int dst_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>xeqScheduleGhostFillNoCoarse</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a30f9669d0d2e51366e08f12a52651e41</anchor>
      <arglist>(int dst_idx, int dst_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>xeqScheduleDataSynch</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ace2d3191f27c726a80a425d998b7a685</anchor>
      <arglist>(int dst_idx, int dst_ln)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt;</type>
      <name>getLevelSAMRAIVectorReal</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>adc73f32eef15519d79d959ea5d9ffc11</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;vec, int level_num) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CCPoissonSolverManager</name>
    <filename>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</filename>
    <member kind="typedef">
      <type>SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;(*</type>
      <name>SolverMaker</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>a7985263917ff1478b5fde862477aea6e</anchor>
      <arglist>)(const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;</type>
      <name>allocateSolver</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>ac61170b4566191465f45a4a0b9ee658c</anchor>
      <arglist>(const std::string &amp;solver_type, const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;</type>
      <name>allocateSolver</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>a3af88e530dbca135481ad0248a54c082</anchor>
      <arglist>(const std::string &amp;solver_type, const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix, const std::string &amp;precond_type, const std::string &amp;precond_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; precond_input_db, const std::string &amp;precond_default_options_prefix) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerSolverFactoryFunction</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>ab2d1fea2517e4407a499a7de0f578f05</anchor>
      <arglist>(const std::string &amp;solver_type, SolverMaker solver_maker)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static CCPoissonSolverManager *</type>
      <name>getManager</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>a6d87c9ab5d200eaf6498b201a85d1dfc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>freeManager</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>ab9bd656fd71da798514f54f1343edd25</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>UNDEFINED</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>aebefb3d01538533931b4e47fb642cf00</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT_KRYLOV_SOLVER</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>a2d6931a66b2e99a525d8d0e5b8adb6ff</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT_FAC_PRECONDITIONER</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>afd08e2e882c8983dcf3c9fc5e2bd7786</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT_LEVEL_SOLVER</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>a35655d823a4f42b277a6e208f5585576</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>CCPoissonSolverManager</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>a1c91ecb22d2f4cbcaec6ee47e3a90e2f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>~CCPoissonSolverManager</name>
      <anchorfile>class_i_b_t_k_1_1_c_c_poisson_solver_manager.html</anchorfile>
      <anchor>a8a4a45788aeaf21a80b8e003bf9181e2</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>IBTK::CellIndexFortranOrder</name>
    <filename>struct_i_b_t_k_1_1_cell_index_fortran_order.html</filename>
    <base>binary_function&lt; SAMRAI::pdat::CellIndex&lt; NDIM &gt;, SAMRAI::pdat::CellIndex&lt; NDIM &gt;, bool &gt;</base>
  </compound>
  <compound kind="class">
    <name>IBTK::CellNoCornersFillPattern</name>
    <filename>class_i_b_t_k_1_1_cell_no_corners_fill_pattern.html</filename>
    <base>VariableFillPattern&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>CellNoCornersFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_cell_no_corners_fill_pattern.html</anchorfile>
      <anchor>a7abc5c689b3eb63f1dd9fea0763141fa</anchor>
      <arglist>(int stencil_width, bool include_dst_patch_box, bool include_edges_on_dst_level, bool include_edges_on_src_level)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CellNoCornersFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_cell_no_corners_fill_pattern.html</anchorfile>
      <anchor>adc4ee533cd91478e409a194adc748d35</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BoxOverlap&lt; NDIM &gt; &gt;</type>
      <name>calculateOverlap</name>
      <anchorfile>class_i_b_t_k_1_1_cell_no_corners_fill_pattern.html</anchorfile>
      <anchor>a2e122df1d835aef77dd07cdb0f735641</anchor>
      <arglist>(const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;dst_geometry, const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;src_geometry, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;dst_patch_box, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;src_mask, bool overwrite_interior, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BoxOverlap&lt; NDIM &gt; &gt;</type>
      <name>calculateOverlapOnLevel</name>
      <anchorfile>class_i_b_t_k_1_1_cell_no_corners_fill_pattern.html</anchorfile>
      <anchor>a76c480d64541f15153be63bf20dd47d4</anchor>
      <arglist>(const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;dst_geometry, const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;src_geometry, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;dst_patch_box, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;src_mask, bool overwrite_interior, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset, int dst_level_num, int src_level_num) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetPatchLevelNumber</name>
      <anchorfile>class_i_b_t_k_1_1_cell_no_corners_fill_pattern.html</anchorfile>
      <anchor>a1f7d2d7dcbbf29656ef65eb9d2856a48</anchor>
      <arglist>(int level_num)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_cell_no_corners_fill_pattern.html</anchorfile>
      <anchor>a956fc5e7ecfaef182b6db46ee7020406</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getPatternName</name>
      <anchorfile>class_i_b_t_k_1_1_cell_no_corners_fill_pattern.html</anchorfile>
      <anchor>a64ccbd10e59d2a51df8a7b72b551989f</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CoarseFineBoundaryRefinePatchStrategy</name>
    <filename>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</filename>
    <base>RefinePatchStrategy&lt; NDIM &gt;</base>
    <member kind="function" virtualness="pure">
      <type>void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>a338e523d31a06706a06bed2ab2bcc2b5</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getRefineOpStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>aa63f000dfc673b08a1a04d476e209980</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>void</type>
      <name>preprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>ac855ed9d5625558f4c58a817492a7ed1</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>void</type>
      <name>postprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>a97be8e50604fc1b14185e12f6f3894c8</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setConsistentInterpolationScheme</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>adf47900f5610ce2b943d62f8589b2e40</anchor>
      <arglist>(bool consistent_type_2_bdry)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setPatchDataIndex</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>a93c99147655de99c47dbf59976f87981</anchor>
      <arglist>(int patch_data_index)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>a51177e7e11d65b5997c2ab1fd69ec9cc</anchor>
      <arglist>(const std::set&lt; int &gt; &amp;patch_data_indices)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setPatchDataIndices</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>af78e1121fe2d40c92b743ae114a0e836</anchor>
      <arglist>(const SAMRAI::hier::ComponentSelector &amp;patch_data_indices)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>aa316651963c517e64554691063cf16fd</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>clearPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>a05a50266dfa79440c8aab6827f2d98fa</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>computeNormalExtension</name>
      <anchorfile>class_i_b_t_k_1_1_coarse_fine_boundary_refine_patch_strategy.html</anchorfile>
      <anchor>a2593302c8125bc1b79d1f2e9dd3a5caf</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CoarsenPatchStrategySet</name>
    <filename>class_i_b_t_k_1_1_coarsen_patch_strategy_set.html</filename>
    <base>CoarsenPatchStrategy&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>CoarsenPatchStrategySet</name>
      <anchorfile>class_i_b_t_k_1_1_coarsen_patch_strategy_set.html</anchorfile>
      <anchor>a350aa09f20a05d6600c57cd35a734f4c</anchor>
      <arglist>(InputIterator first, InputIterator last, bool managed=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CoarsenPatchStrategySet</name>
      <anchorfile>class_i_b_t_k_1_1_coarsen_patch_strategy_set.html</anchorfile>
      <anchor>af0bac321cbf79b8b781d0ee18890be64</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getCoarsenOpStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_coarsen_patch_strategy_set.html</anchorfile>
      <anchor>a4da8b2848a8cd8de43f4dcfb52d79f01</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessCoarsen</name>
      <anchorfile>class_i_b_t_k_1_1_coarsen_patch_strategy_set.html</anchorfile>
      <anchor>a7c1e63a6f29d429bc948a66595d9537a</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;coarse_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessCoarsen</name>
      <anchorfile>class_i_b_t_k_1_1_coarsen_patch_strategy_set.html</anchorfile>
      <anchor>abb5a08830af8764cc72b0f2a5ce6ea2a</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;coarse_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CopyToRootSchedule</name>
    <filename>class_i_b_t_k_1_1_copy_to_root_schedule.html</filename>
    <member kind="function">
      <type></type>
      <name>CopyToRootSchedule</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_schedule.html</anchorfile>
      <anchor>acfde55bc0ee25cf358cac414ccea5c72</anchor>
      <arglist>(int root_proc, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level, int src_patch_data_idx)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CopyToRootSchedule</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_schedule.html</anchorfile>
      <anchor>a3416872f891c8b64686624875053178b</anchor>
      <arglist>(int root_proc, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level, const std::vector&lt; int &gt; &amp;src_patch_data_idxs)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CopyToRootSchedule</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_schedule.html</anchorfile>
      <anchor>a28e08da57141966af076a525957e8f51</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>communicate</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_schedule.html</anchorfile>
      <anchor>a8fba46c3bb1a6e844aa572de446c0136</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchData&lt; NDIM &gt; &gt; &gt; &amp;</type>
      <name>getRootPatchData</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_schedule.html</anchorfile>
      <anchor>acd8e988e66dba9cc8f1369552c3da494</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::CopyToRootTransaction</name>
    <filename>class_i_b_t_k_1_1_copy_to_root_transaction.html</filename>
    <base>SAMRAI::tbox::Transaction</base>
    <member kind="function">
      <type></type>
      <name>CopyToRootTransaction</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>a1128518560195d0544a6403151ad3ef0</anchor>
      <arglist>(int src_proc, int dst_proc, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level, int src_patch_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchData&lt; NDIM &gt; &gt; dst_patch_data)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CopyToRootTransaction</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>a5ad2e080be30c7101a01cca85d28e400</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchData&lt; NDIM &gt; &gt;</type>
      <name>getRootPatchData</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>a643f7d9194baed4409e54624683ab5a1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>canEstimateIncomingMessageSize</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>acad549672fa54f4172d16739a9b7b9b3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>computeIncomingMessageSize</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>ac761d6f9ae03cb1c7a48320a4c52783b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>computeOutgoingMessageSize</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>a37a13f062298752bdbcf8c1b35c69da6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getSourceProcessor</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>a7454d1896ac996dd6a3b778319cd338c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getDestinationProcessor</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>a311b0c1fd3451f6874ba9a2e461bb9b5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>a5d778a5730ae0bb24501174b695f14c6</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unpackStream</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>af749928697e475c5925f17c395b81460</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>copyLocalData</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>a994e2941b4fdac1aa89813889e28ccc6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>printClassData</name>
      <anchorfile>class_i_b_t_k_1_1_copy_to_root_transaction.html</anchorfile>
      <anchor>aa6df04fbefb0bb45ccb27195b16ef8ea</anchor>
      <arglist>(std::ostream &amp;stream) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::DebuggingUtilities</name>
    <filename>class_i_b_t_k_1_1_debugging_utilities.html</filename>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>checkCellDataForNaNs</name>
      <anchorfile>class_i_b_t_k_1_1_debugging_utilities.html</anchorfile>
      <anchor>a9958c6427c2f24fcb60e830af7708d6d</anchor>
      <arglist>(int patch_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, bool interior_only=true, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>checkFaceDataForNaNs</name>
      <anchorfile>class_i_b_t_k_1_1_debugging_utilities.html</anchorfile>
      <anchor>ade61235cfb06f9425f8e11875355ed37</anchor>
      <arglist>(int patch_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, bool interior_only=true, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>checkNodeDataForNaNs</name>
      <anchorfile>class_i_b_t_k_1_1_debugging_utilities.html</anchorfile>
      <anchor>a017c406aa867147ea4e982beca923f70</anchor>
      <arglist>(int patch_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, bool interior_only=true, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>checkSideDataForNaNs</name>
      <anchorfile>class_i_b_t_k_1_1_debugging_utilities.html</anchorfile>
      <anchor>afdcc2502273f7d2ee9487c7af4123ff1</anchor>
      <arglist>(int patch_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, bool interior_only=true, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>saveCellData</name>
      <anchorfile>class_i_b_t_k_1_1_debugging_utilities.html</anchorfile>
      <anchor>ad34b913fb71daec0c4b9c81e6a18784f</anchor>
      <arglist>(int patch_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, const std::string &amp;filename, const std::string &amp;dirname)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>saveFaceData</name>
      <anchorfile>class_i_b_t_k_1_1_debugging_utilities.html</anchorfile>
      <anchor>a64afd2e85c4bb296a3f808586ca64507</anchor>
      <arglist>(int patch_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, const std::string &amp;filename, const std::string &amp;dirname)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>saveNodeData</name>
      <anchorfile>class_i_b_t_k_1_1_debugging_utilities.html</anchorfile>
      <anchor>a53f621cd5fe9aac6f866692df9b40b6c</anchor>
      <arglist>(int patch_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, const std::string &amp;filename, const std::string &amp;dirname)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>saveSideData</name>
      <anchorfile>class_i_b_t_k_1_1_debugging_utilities.html</anchorfile>
      <anchor>a117f384942d8e37a3a4ee87c8b986e38</anchor>
      <arglist>(int patch_data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, const std::string &amp;filename, const std::string &amp;dirname)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>saveLagrangianData</name>
      <anchorfile>class_i_b_t_k_1_1_debugging_utilities.html</anchorfile>
      <anchor>a46fa27f48834e7504c0e938ffe77ed6e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LData &gt; lag_data, bool save_ghost_nodes, const std::string &amp;filename, const std::string &amp;dirname)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>IBTK::DofObjectComp</name>
    <filename>struct_i_b_t_k_1_1_dof_object_comp.html</filename>
    <base>binary_function&lt; const libMesh::DofObject *const, const libMesh::DofObject *const, bool &gt;</base>
  </compound>
  <compound kind="class">
    <name>IBTK::EdgeDataSynchronization</name>
    <filename>class_i_b_t_k_1_1_edge_data_synchronization.html</filename>
    <class kind="class">IBTK::EdgeDataSynchronization::SynchronizationTransactionComponent</class>
    <member kind="function">
      <type></type>
      <name>EdgeDataSynchronization</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization.html</anchorfile>
      <anchor>abe3fbf4786188659e1a7f311cdfda26a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~EdgeDataSynchronization</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization.html</anchorfile>
      <anchor>a6f6a3da6575f3d613c1c389a35414342</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization.html</anchorfile>
      <anchor>a526bc2a734c0a18341302692dccd660d</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;transaction_comp, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization.html</anchorfile>
      <anchor>a416778340c4caad32959f15bd16124c5</anchor>
      <arglist>(const std::vector&lt; SynchronizationTransactionComponent &gt; &amp;transaction_comps, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization.html</anchorfile>
      <anchor>aac264a73d7bf1f654f3158e795585d98</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;transaction_comps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTransactionComponents</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization.html</anchorfile>
      <anchor>ac3e14b54d1f50991ecdebefebbff62a3</anchor>
      <arglist>(const std::vector&lt; SynchronizationTransactionComponent &gt; &amp;transaction_comps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization.html</anchorfile>
      <anchor>a99560d6e56507048ec61e1a68d067903</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>synchronizeData</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization.html</anchorfile>
      <anchor>a6417b7870610be3be5d90c8beecfe131</anchor>
      <arglist>(double fill_time)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::EdgeDataSynchronization::SynchronizationTransactionComponent</name>
    <filename>class_i_b_t_k_1_1_edge_data_synchronization_1_1_synchronization_transaction_component.html</filename>
    <member kind="function">
      <type></type>
      <name>SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>aa874b4a203e5357f153fb3390cdf6ca3</anchor>
      <arglist>(int data_idx=-1, const std::string &amp;coarsen_op_name=&quot;NONE&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>acaf626b83c4d1b6b04c2cb7097957451</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;from)</arglist>
    </member>
    <member kind="function">
      <type>SynchronizationTransactionComponent &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>a2f48b64de40ca6c82d1cb6f213f8a29a</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;that)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_edge_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>aacd201a960fb33145bba5d62181522dc</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::EdgeSynchCopyFillPattern</name>
    <filename>class_i_b_t_k_1_1_edge_synch_copy_fill_pattern.html</filename>
    <base>VariableFillPattern&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>EdgeSynchCopyFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_edge_synch_copy_fill_pattern.html</anchorfile>
      <anchor>accbf9d4302d02fb745dd9ea7f5b7d7cc</anchor>
      <arglist>(unsigned int axis)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~EdgeSynchCopyFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_edge_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a801c60a80e74abd12efc5b20b70fdb7a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BoxOverlap&lt; NDIM &gt; &gt;</type>
      <name>calculateOverlap</name>
      <anchorfile>class_i_b_t_k_1_1_edge_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a169b5d3d7b3401393400326953a7496e</anchor>
      <arglist>(const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;dst_geometry, const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;src_geometry, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;dst_patch_box, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;src_mask, bool overwrite_interior, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_edge_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a03cd7300427e2bca913b85aaf0758a2d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getPatternName</name>
      <anchorfile>class_i_b_t_k_1_1_edge_synch_copy_fill_pattern.html</anchorfile>
      <anchor>ab6926c48ac0312171c4e0de0679012bc</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::ExtendedRobinBcCoefStrategy</name>
    <filename>class_i_b_t_k_1_1_extended_robin_bc_coef_strategy.html</filename>
    <base>RobinBcCoefStrategy&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>ExtendedRobinBcCoefStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_extended_robin_bc_coef_strategy.html</anchorfile>
      <anchor>a9569e51eb89192b505a05bcd550d2afb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ExtendedRobinBcCoefStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_extended_robin_bc_coef_strategy.html</anchorfile>
      <anchor>aaf4f35c87d8bace752b1fe5ca2f3a95e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setTargetPatchDataIndex</name>
      <anchorfile>class_i_b_t_k_1_1_extended_robin_bc_coef_strategy.html</anchorfile>
      <anchor>a4c03093ecfa83604a1083dfc03c8a443</anchor>
      <arglist>(int target_data_idx)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>clearTargetPatchDataIndex</name>
      <anchorfile>class_i_b_t_k_1_1_extended_robin_bc_coef_strategy.html</anchorfile>
      <anchor>a0a2445f08a61823fde61d86e1708a7ee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_extended_robin_bc_coef_strategy.html</anchorfile>
      <anchor>a9f8ae787f9c0600c32470ea78c16c2dc</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::FaceDataSynchronization</name>
    <filename>class_i_b_t_k_1_1_face_data_synchronization.html</filename>
    <class kind="class">IBTK::FaceDataSynchronization::SynchronizationTransactionComponent</class>
    <member kind="function">
      <type></type>
      <name>FaceDataSynchronization</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization.html</anchorfile>
      <anchor>a38793dd36d7890a721122edd81114f91</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~FaceDataSynchronization</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization.html</anchorfile>
      <anchor>aba03337dcba3df2d910f959920164261</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization.html</anchorfile>
      <anchor>a05a1d862cbcf5984ca78077fc22422dc</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;transaction_comp, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization.html</anchorfile>
      <anchor>a59d142c739f94edc5c487e23cf180f02</anchor>
      <arglist>(const std::vector&lt; SynchronizationTransactionComponent &gt; &amp;transaction_comps, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization.html</anchorfile>
      <anchor>a43e798624ac28c8d8f24102c64dc0a7f</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;transaction_comps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTransactionComponents</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization.html</anchorfile>
      <anchor>ad88827e72ad94a7dfbab31cf7e86f9a5</anchor>
      <arglist>(const std::vector&lt; SynchronizationTransactionComponent &gt; &amp;transaction_comps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization.html</anchorfile>
      <anchor>a3a4cc04666817ea280b526ddcd5de98b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>synchronizeData</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization.html</anchorfile>
      <anchor>a2d354c4f41e211509efc59975eb424a8</anchor>
      <arglist>(double fill_time)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::FaceDataSynchronization::SynchronizationTransactionComponent</name>
    <filename>class_i_b_t_k_1_1_face_data_synchronization_1_1_synchronization_transaction_component.html</filename>
    <member kind="function">
      <type></type>
      <name>SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>aaab23d4fa056d92a2f5b7af57d37a7dc</anchor>
      <arglist>(int data_idx=-1, const std::string &amp;coarsen_op_name=&quot;NONE&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>adb2f586994cdf5d5b92a16a21fd4acd0</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;from)</arglist>
    </member>
    <member kind="function">
      <type>SynchronizationTransactionComponent &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>a39f631f460a5bf16efc66930ce46173c</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;that)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_face_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>ababb30b66cff33c29298cafb73a6f37e</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::FaceSynchCopyFillPattern</name>
    <filename>class_i_b_t_k_1_1_face_synch_copy_fill_pattern.html</filename>
    <base>VariableFillPattern&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>FaceSynchCopyFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_face_synch_copy_fill_pattern.html</anchorfile>
      <anchor>ad3669abcc4a39877cf0eda41767d7d51</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~FaceSynchCopyFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_face_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a507e7d2697779446b875cb0859b3b3b1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BoxOverlap&lt; NDIM &gt; &gt;</type>
      <name>calculateOverlap</name>
      <anchorfile>class_i_b_t_k_1_1_face_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a0856edc9fa86cee52058576b2e15f7af</anchor>
      <arglist>(const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;dst_geometry, const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;src_geometry, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;dst_patch_box, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;src_mask, bool overwrite_interior, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_face_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a1840d17772b72783b9afd8754070c7b6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getPatternName</name>
      <anchorfile>class_i_b_t_k_1_1_face_synch_copy_fill_pattern.html</anchorfile>
      <anchor>af2b7f6d1567c43af4e081be9402ecd2a</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::FACPreconditioner</name>
    <filename>class_i_b_t_k_1_1_f_a_c_preconditioner.html</filename>
    <base>IBTK::LinearSolver</base>
    <member kind="function">
      <type></type>
      <name>FACPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a424005fe128542abb73f05037458c085</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; FACPreconditionerStrategy &gt; fac_strategy, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~FACPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a2e605be84c5ce70049fec358ec8870a0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a6bdd459f5d29f33c9792579cdae33a1e</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSolutionTime</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a8e1117b9a0d2d36cbb44a1d177148594</anchor>
      <arglist>(double solution_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a953b7b1d13ff0cd0836f3e76f9b46070</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>aed7649ecfc99a82bbd0035b67d2d110a</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a5ce4e5832843371020202ecee696769a</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>ada0e65b2016918daed7861cbc24bf7c3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setInitialGuessNonzero</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>aff46f7f8d8c2f4c22ec97133ffcee524</anchor>
      <arglist>(bool initial_guess_nonzero=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMaxIterations</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a16e89c860f8010fe1510a9cd8f374ee8</anchor>
      <arglist>(int max_iterations)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMGCycleType</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a21abb475a19ae8b44cbac572caa3154a</anchor>
      <arglist>(MGCycleType cycle_type)</arglist>
    </member>
    <member kind="function">
      <type>MGCycleType</type>
      <name>getMGCycleType</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a05138be367c0013a977edeee0bb11d15</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setNumPreSmoothingSweeps</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a01f25b334c9e91120606e4d86ee666b4</anchor>
      <arglist>(int num_pre_sweeps)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumPreSmoothingSweeps</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a3965149020d6931cad3719ce578c51b7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setNumPostSmoothingSweeps</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>af66858d37fc635c791318a2692b22870</anchor>
      <arglist>(int num_post_sweeps)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumPostSmoothingSweeps</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner.html</anchorfile>
      <anchor>a5f033a916d7248a22051f79c26a2d0d0</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::FACPreconditionerStrategy</name>
    <filename>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>restrictResidual</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a168f159bc2cad7c3341221908fbb4216</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;source, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;dest, int dest_level_num)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>prolongError</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a107c5d2ca956ebd015dee02de6cc7f92</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;source, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;dest, int dest_level_num)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>prolongErrorAndCorrect</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ab74de4f4dca3a0442208645465139e7a</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;source, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;dest, int dest_level_num)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>smoothError</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ab0b305ed685e65b6e1e9588e68a63d72</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;error, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, int level_num, int num_sweeps, bool performing_pre_sweeps, bool performing_post_sweeps)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>solveCoarsestLevel</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ad8bbbc83a58aa680b19cb84605032a00</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;error, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, int coarsest_level_num)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>computeResidual</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a0e3578c0f7cad8a598c676bf47fde61b</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs, int coarsest_level_num, int finest_level_num)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a7ef0ff73722786c412ef26e56d46831c</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a910430e83b5309d2fb8051c497772f22</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>allocateScratchData</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>aa9b16d759b9cbe31411d2e75dfe9a54d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>deallocateScratchData</name>
      <anchorfile>class_i_b_t_k_1_1_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>afa2dbd0a0a59eef35604173b4ef26750</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::FEDataManager</name>
    <filename>class_i_b_t_k_1_1_f_e_data_manager.html</filename>
    <base>SAMRAI::tbox::Serializable</base>
    <base>StandardTagAndInitStrategy&lt; NDIM &gt;</base>
    <class kind="struct">IBTK::FEDataManager::InterpSpec</class>
    <class kind="struct">IBTK::FEDataManager::SpreadSpec</class>
    <member kind="function">
      <type>void</type>
      <name>registerLoadBalancer</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>ab1a0a70fa3d4a343a8bdb63b7f89b9fb</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::LoadBalancer&lt; NDIM &gt; &gt; load_balancer, int workload_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setEquationSystems</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a7ab2051403b2a6d0038cb5792c0017e1</anchor>
      <arglist>(libMesh::EquationSystems *equation_systems, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>libMesh::EquationSystems *</type>
      <name>getEquationSystems</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a75fdc4552da15952615ecfadefb19fa1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getLevelNumber</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a03eeee1f8a83928b7454a254c13fcce4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getGhostCellWidth</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a54a0c4134408f7d0c2027f70539c7d15</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const InterpSpec &amp;</type>
      <name>getDefaultInterpSpec</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a9f00f57c1f18f06ed9c4cbe5def7caf9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const SpreadSpec &amp;</type>
      <name>getDefaultSpreadSpec</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a7c6d96b29c80d107f69adf1fe0df411a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; std::vector&lt; libMesh::Elem * &gt; &gt; &amp;</type>
      <name>getActivePatchElementMap</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>aeb37f9460e1b45d7d8e454cd344f6577</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinitElementMappings</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a3d232e10e554539a93a127e02a48dc21</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>libMesh::NumericVector&lt; double &gt; *</type>
      <name>getSolutionVector</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>aa6bcac630c5d01de43874c99f716ccf6</anchor>
      <arglist>(const std::string &amp;system_name) const </arglist>
    </member>
    <member kind="function">
      <type>libMesh::NumericVector&lt; double &gt; *</type>
      <name>buildGhostedSolutionVector</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a28ee72e7fdd4007414c9f7f44afa2806</anchor>
      <arglist>(const std::string &amp;system_name, bool localize_data=true)</arglist>
    </member>
    <member kind="function">
      <type>libMesh::NumericVector&lt; double &gt; *</type>
      <name>getCoordsVector</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a69c2939be6577b4f4061b7fa196e37b8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>libMesh::NumericVector&lt; double &gt; *</type>
      <name>buildGhostedCoordsVector</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>aea3e740109fb669d08420a059720028b</anchor>
      <arglist>(bool localize_data=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>af4cf57b4aa661f3f29c343f979ff264c</anchor>
      <arglist>(int f_data_idx, libMesh::NumericVector&lt; double &gt; &amp;F, libMesh::NumericVector&lt; double &gt; &amp;X, const std::string &amp;system_name, RobinPhysBdryPatchStrategy *f_phys_bdry_op, double fill_data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>aadcfbaf4a4ce2dc7c66931e89c26ada1</anchor>
      <arglist>(int f_data_idx, libMesh::NumericVector&lt; double &gt; &amp;F, libMesh::NumericVector&lt; double &gt; &amp;X, const std::string &amp;system_name, const SpreadSpec &amp;spread_spec, RobinPhysBdryPatchStrategy *f_phys_bdry_op, double fill_data_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>prolongData</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a215cd4386b8467e60b183573cdb1fdd8</anchor>
      <arglist>(int f_data_idx, libMesh::NumericVector&lt; double &gt; &amp;F, libMesh::NumericVector&lt; double &gt; &amp;X, const std::string &amp;system_name, bool is_density=true, bool accumulate_on_grid=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>ac869835ef393917f0d5202d4bb644c6d</anchor>
      <arglist>(int f_data_idx, libMesh::NumericVector&lt; double &gt; &amp;F, libMesh::NumericVector&lt; double &gt; &amp;X, const std::string &amp;system_name, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_refine_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>ab3e9302c3fc4a37d3da6f4f8137b42fe</anchor>
      <arglist>(int f_data_idx, libMesh::NumericVector&lt; double &gt; &amp;F, libMesh::NumericVector&lt; double &gt; &amp;X, const std::string &amp;system_name, const InterpSpec &amp;interp_spec, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_refine_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>restrictData</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>ab2bdfdf74517773f8ef216c3a9a7180e</anchor>
      <arglist>(int f_data_idx, libMesh::NumericVector&lt; double &gt; &amp;F, libMesh::NumericVector&lt; double &gt; &amp;X, const std::string &amp;system_name, bool use_consistent_mass_matrix=true)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; libMesh::LinearSolver&lt; double &gt; *, libMesh::SparseMatrix&lt; double &gt; * &gt;</type>
      <name>buildL2ProjectionSolver</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>abb0a82be8dd81dd70af8d1e89ca2d1dd</anchor>
      <arglist>(const std::string &amp;system_name, libMeshEnums::QuadratureType quad_type=QGAUSS, libMeshEnums::Order quad_order=FIFTH)</arglist>
    </member>
    <member kind="function">
      <type>libMesh::NumericVector&lt; double &gt; *</type>
      <name>buildDiagonalL2MassMatrix</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>ac8451f2797a06ffeb2f2d4c49bb3bb65</anchor>
      <arglist>(const std::string &amp;system_name)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>computeL2Projection</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a4f8a80bbbe3f0454ce00bf0da8c61218</anchor>
      <arglist>(libMesh::NumericVector&lt; double &gt; &amp;U, libMesh::NumericVector&lt; double &gt; &amp;F, const std::string &amp;system_name, bool consistent_mass_matrix=true, libMeshEnums::QuadratureType quad_type=QGAUSS, libMeshEnums::Order quad_order=FIFTH, double tol=1.0e-6, unsigned int max_its=100)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>updateWorkloadEstimates</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>aa356523e0dbbf1f259cb678b09d1fee9</anchor>
      <arglist>(int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a2afd491723bc5b1daad8f220e4f9d0b1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt;(NULL), bool allocate_data=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetHierarchyConfiguration</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a79b27d1aff937f2ce5e8aa84424187a0</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_ln, int finest_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyGradientDetector</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a6e2532c268ca81bf0c31ee95298d22c1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a24178a8d43bf01b2394eeb243e44df9c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a7cd0c81a11489e5d334cafa744253e70</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt;</type>
      <name>getPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>ae3f02bf198101e2b4f716bb2537b8df0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchLevels</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>ad642f1a5b469f6fa128a97f843fb91cd</anchor>
      <arglist>(int coarsest_ln, int finest_ln)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; int, int &gt;</type>
      <name>getPatchLevels</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>ac498f8acc7c2c93532f7dfc51933e1ae</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" static="yes">
      <type>static FEDataManager *</type>
      <name>getManager</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a0f1758b25d2c421aee5968036031189d</anchor>
      <arglist>(const std::string &amp;name, const InterpSpec &amp;default_interp_spec, const SpreadSpec &amp;default_spread_spec, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;min_ghost_width=SAMRAI::hier::IntVector&lt; NDIM &gt;(0), bool register_for_restart=true)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>freeAllManagers</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a88429e0a8a395f3552efeb71de52a01a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>updateInterpQuadratureRule</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a520171a0d212575eadea75914246b547</anchor>
      <arglist>(libMesh::AutoPtr&lt; libMesh::QBase &gt; &amp;qrule, const InterpSpec &amp;spec, libMesh::Elem *elem, const boost::multi_array&lt; double, 2 &gt; &amp;X_node, double dx_min)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>updateSpreadQuadratureRule</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>ad29c8e988c8c995a6b38d1e992b4b561</anchor>
      <arglist>(libMesh::AutoPtr&lt; libMesh::QBase &gt; &amp;qrule, const SpreadSpec &amp;spec, libMesh::Elem *elem, const boost::multi_array&lt; double, 2 &gt; &amp;X_node, double dx_min)</arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>COORDINATES_SYSTEM_NAME</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a1a074b6d7fdbc4e48735e970b33100f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const short int</type>
      <name>ZERO_DISPLACEMENT_X_BDRY_ID</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a3e76a91c964c5e5ed8f631a261b8ed3d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>FEDataManager</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>a8e38c56f46d2f07c090187f99257e29a</anchor>
      <arglist>(const std::string &amp;object_name, const InterpSpec &amp;default_interp_spec, const SpreadSpec &amp;default_spread_spec, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>~FEDataManager</name>
      <anchorfile>class_i_b_t_k_1_1_f_e_data_manager.html</anchorfile>
      <anchor>afefddb2a11c0a9c77fdbcdf6f6be6896</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>IBTK::FEDataManager::InterpSpec</name>
    <filename>struct_i_b_t_k_1_1_f_e_data_manager_1_1_interp_spec.html</filename>
  </compound>
  <compound kind="struct">
    <name>IBTK::FEDataManager::SpreadSpec</name>
    <filename>struct_i_b_t_k_1_1_f_e_data_manager_1_1_spread_spec.html</filename>
  </compound>
  <compound kind="class">
    <name>IBTK::FixedSizedStream</name>
    <filename>class_i_b_t_k_1_1_fixed_sized_stream.html</filename>
    <base>SAMRAI::tbox::AbstractStream</base>
    <member kind="function">
      <type></type>
      <name>FixedSizedStream</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a45dc421969d89e561830d1168ea53656</anchor>
      <arglist>(int bytes)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FixedSizedStream</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>ae95a720acf3c6524891da47c772c52ee</anchor>
      <arglist>(const void *buffer, int bytes)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~FixedSizedStream</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>aa8467cda32148cb4f93d9a41ea19aaee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>getBufferStart</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a99e36e606de6ef10d291c118bf34bc94</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const void *</type>
      <name>getBufferStart</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>abf70e56327252aedcc446d75c871ebed</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getCurrentSize</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>ae032bdc682f88e8bfdc33b1e8ecdd916</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getCurrentIndex</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>ac87fe71a1727fdba1b9ab9f247ca80c3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCurrentIndex</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>aa54e5ec695adb461e8863b746d17dc39</anchor>
      <arglist>(int index)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetIndex</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a1dc641a90eb7aacd522798760942bd9d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>affba8bfffe801be08ebdb69dc2a96a6f</anchor>
      <arglist>(const bool &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>afa47be0f5f6c3de34f9843213b4ab6a5</anchor>
      <arglist>(bool &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>ae5b865c00b1d3b03c6d1c29162bad19b</anchor>
      <arglist>(const bool *data, int n=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unpack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>aacf4eae0a404b68196ea0b77cdc22a03</anchor>
      <arglist>(bool *data, int n=1)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a56b282d906c1452887218f8c42b6c8ba</anchor>
      <arglist>(const char &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>aabceb0dabf61c45db6ebf9abbba0c773</anchor>
      <arglist>(char &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a9db778c6b6f7bb68c9fcbbb51b0647b8</anchor>
      <arglist>(const char *data, int n=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unpack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a1ab67eaffd9cd938cf71f1eacb9af31d</anchor>
      <arglist>(char *data, int n=1)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a00ace1ed47d9051010fcc1e880f24cb0</anchor>
      <arglist>(const dcomplex &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a32446e01dbe760cd87f10169b2510a03</anchor>
      <arglist>(dcomplex &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>abb22814a05b34b63ddcc14039156c43c</anchor>
      <arglist>(const dcomplex *data, int n=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unpack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a59375a319336ae64884b26385470b552</anchor>
      <arglist>(dcomplex *data, int n=1)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>af429e729bd3d15ba03605a64a2b6cb1c</anchor>
      <arglist>(const double &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a44fbe58be1d40d8385877487e52535bd</anchor>
      <arglist>(double &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a0d476856e2ec16e31bea1cd95d69762d</anchor>
      <arglist>(const double *data, int n=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unpack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>acfb1e210ca0fbccaab75c3c46ddbf1b4</anchor>
      <arglist>(double *data, int n=1)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a73dcb8a16df6dcd8ccc41372a8388787</anchor>
      <arglist>(const float &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>acbf7045798972f64c079317339980467</anchor>
      <arglist>(float &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a97fefc28c73cfa992bc6754d3da14751</anchor>
      <arglist>(const float *data, int n=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unpack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>aa1869e675480e4a804fb8e34fc702123</anchor>
      <arglist>(float *data, int n=1)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>aa32239f86419166b3b096e3f58af4437</anchor>
      <arglist>(const int &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::AbstractStream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a9a809baee35b9d0f9457aa9e97f8e5b3</anchor>
      <arglist>(int &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>a1dad471e72289ecbcee6f018cad2855b</anchor>
      <arglist>(const int *data, int n=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unpack</name>
      <anchorfile>class_i_b_t_k_1_1_fixed_sized_stream.html</anchorfile>
      <anchor>aab2599db90de8db59e89cfc3e1198c9c</anchor>
      <arglist>(int *data, int n=1)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::GeneralOperator</name>
    <filename>class_i_b_t_k_1_1_general_operator.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>apply</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a29d5dfa5f65e28ff7f52866a6162e6eb</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>a2c4852d470d984e34bd15d7b0dbb5418</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_general_operator.html</anchorfile>
      <anchor>af1dc37c17a2a2b12e67ad88f0a5d6a6b</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::GeneralSolver</name>
    <filename>class_i_b_t_k_1_1_general_solver.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a8c006d4a74c8860ef66b095ee93a51e6</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>ac88664e72b915b797564e49e5520edb0</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>a0d0d26dc3df7083c533663bef216cd0f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>printClassData</name>
      <anchorfile>class_i_b_t_k_1_1_general_solver.html</anchorfile>
      <anchor>aa775ba22b544df90645d15f193bcf2b5</anchor>
      <arglist>(std::ostream &amp;stream)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::HierarchyGhostCellInterpolation</name>
    <filename>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation.html</filename>
    <class kind="class">IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent</class>
    <member kind="function">
      <type></type>
      <name>HierarchyGhostCellInterpolation</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation.html</anchorfile>
      <anchor>a602137b35713dd6aa3da83a518aee748</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~HierarchyGhostCellInterpolation</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation.html</anchorfile>
      <anchor>accde30040f037ac1c29c07690f56a3d5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation.html</anchorfile>
      <anchor>a174ad5213c2fb573a0b537a84b6319c0</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation.html</anchorfile>
      <anchor>a19ec30ce9549ab86009a526b79f55066</anchor>
      <arglist>(InterpolationTransactionComponent transaction_comp, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation.html</anchorfile>
      <anchor>a23aac6e6a790824d1c1b5be0ad9774c4</anchor>
      <arglist>(const std::vector&lt; InterpolationTransactionComponent &gt; &amp;transaction_comps, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation.html</anchorfile>
      <anchor>aa0f06871844e77712215d06e2d8b4c3a</anchor>
      <arglist>(const InterpolationTransactionComponent &amp;transaction_comps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTransactionComponents</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation.html</anchorfile>
      <anchor>ae8de868a0933722f56c15b14cf4e9de8</anchor>
      <arglist>(const std::vector&lt; InterpolationTransactionComponent &gt; &amp;transaction_comps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinitializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation.html</anchorfile>
      <anchor>a4e33acdab6fce4bd5b17c6a62ee19625</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation.html</anchorfile>
      <anchor>a5accc85b3939118d459dc848bd29e014</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>fillData</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation.html</anchorfile>
      <anchor>a79c21b336cc4585c26a4368f931df697</anchor>
      <arglist>(double fill_time)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent</name>
    <filename>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation_1_1_interpolation_transaction_component.html</filename>
    <member kind="function">
      <type></type>
      <name>InterpolationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation_1_1_interpolation_transaction_component.html</anchorfile>
      <anchor>a568c66401198a9077bcb14711a137569</anchor>
      <arglist>(int data_idx=-1, const std::string &amp;refine_op_name=&quot;NONE&quot;, bool use_cf_bdry_interpolation=false, const std::string &amp;coarsen_op_name=&quot;NONE&quot;, const std::string &amp;phys_bdry_extrap_type=&quot;NONE&quot;, bool consistent_type_2_bdry=false, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *robin_bc_coef=NULL, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::VariableFillPattern&lt; NDIM &gt; &gt; fill_pattern=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InterpolationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation_1_1_interpolation_transaction_component.html</anchorfile>
      <anchor>a14ab30c5d96a18434c7dcda3ccb60d0c</anchor>
      <arglist>(int data_idx, const std::string &amp;refine_op_name, bool use_cf_bdry_interpolation, const std::string &amp;coarsen_op_name, const std::string &amp;phys_bdry_extrap_type, bool consistent_type_2_bdry, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;robin_bc_coefs, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::VariableFillPattern&lt; NDIM &gt; &gt; fill_pattern=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InterpolationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation_1_1_interpolation_transaction_component.html</anchorfile>
      <anchor>af3a339d325422a667254b5ee898ba820</anchor>
      <arglist>(int dst_data_idx, int src_data_idx, const std::string &amp;refine_op_name, bool use_cf_bdry_interpolation, const std::string &amp;coarsen_op_name, const std::string &amp;phys_bdry_extrap_type, bool consistent_type_2_bdry, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *robin_bc_coef, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::VariableFillPattern&lt; NDIM &gt; &gt; fill_pattern=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InterpolationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation_1_1_interpolation_transaction_component.html</anchorfile>
      <anchor>a8bab7343beb475d7421dd1c5d4a1746b</anchor>
      <arglist>(int dst_data_idx, int src_data_idx, const std::string &amp;refine_op_name, bool use_cf_bdry_interpolation, const std::string &amp;coarsen_op_name, const std::string &amp;phys_bdry_extrap_type, bool consistent_type_2_bdry, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;robin_bc_coefs, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::VariableFillPattern&lt; NDIM &gt; &gt; fill_pattern=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InterpolationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation_1_1_interpolation_transaction_component.html</anchorfile>
      <anchor>a9d456a8d10fc95072777b6e0d7235af0</anchor>
      <arglist>(const InterpolationTransactionComponent &amp;from)</arglist>
    </member>
    <member kind="function">
      <type>InterpolationTransactionComponent &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation_1_1_interpolation_transaction_component.html</anchorfile>
      <anchor>aaaf250082d6424fcaffa22577c06fde4</anchor>
      <arglist>(const InterpolationTransactionComponent &amp;that)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~InterpolationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_ghost_cell_interpolation_1_1_interpolation_transaction_component.html</anchorfile>
      <anchor>ab7eec01f4aa31d06c505fb0aa8a8cde7</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::HierarchyIntegrator</name>
    <filename>class_i_b_t_k_1_1_hierarchy_integrator.html</filename>
    <base>StandardTagAndInitStrategy&lt; NDIM &gt;</base>
    <base>SAMRAI::tbox::Serializable</base>
    <member kind="typedef">
      <type>void(*</type>
      <name>PreprocessIntegrateHierarchyCallbackFcnPtr</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>aebe9971c7d1d3d856391fa32a82eeb4b</anchor>
      <arglist>)(double current_time, double new_time, int num_cycles, void *ctx)</arglist>
    </member>
    <member kind="typedef">
      <type>void(*</type>
      <name>IntegrateHierarchyCallbackFcnPtr</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a838715fb4c2588f60c07bd647bc1bd90</anchor>
      <arglist>)(double current_time, double new_time, int cycle_num, void *ctx)</arglist>
    </member>
    <member kind="typedef">
      <type>void(*</type>
      <name>PostprocessIntegrateHierarchyCallbackFcnPtr</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ab72e9f5f0a34f643e8a329ef437777ad</anchor>
      <arglist>)(double current_time, double new_time, bool skip_synchronize_new_state_data, int num_cycles, void *ctx)</arglist>
    </member>
    <member kind="typedef">
      <type>void(*</type>
      <name>ApplyGradientDetectorCallbackFcnPtr</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a984a22e86b3ca2130c22209f7232a2b6</anchor>
      <arglist>)(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too, void *ctx)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>HierarchyIntegrator</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ab72f87f4a0fae64ef1f323bd015e04f1</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, bool register_for_restart)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~HierarchyIntegrator</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a333eb83be8d87d8eec11a914edb6cf19</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getName</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a86bbd2195d7cb018ce75fdefe1998120</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>initializeHierarchyIntegrator</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a6fe5881b1214ae1796b2c48cbb5943b0</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializePatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a3483a9cd8e7943b6409f184593ad49bd</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>advanceHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ad03e2b9b63c8cc24cd65085a49593d19</anchor>
      <arglist>(double dt)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getMaximumTimeStepSize</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a0d132b788a509557fd1d5c9fe5aa2a4b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>synchronizeHierarchyData</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ab70f9c1d30a3315db71148143297f6b4</anchor>
      <arglist>(VariableContextType ctx_type)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTimeDependentHierarchyData</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a2bf10691b4bf6a42a7206230335b31ab</anchor>
      <arglist>(double new_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetIntegratorToPreadvanceState</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ac779a5411a28b6a42b48f8e23105418c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>regridHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a40f8bd640570540e221243aa29832527</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>atRegridPoint</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a3cc479289b6c0140416d313c8bcf4a44</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getIntegratorTime</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a0322680b9047cdb03a1350fda3296e9b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getStartTime</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a1637678b12e645ab285c6488e63308da</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getEndTime</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>abdd00131c17c5e2b1e321165b689754e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getIntegratorStep</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a63a0f763e51dc6e0e237948c9b1dbc5b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getMaxIntegratorSteps</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>aa19551072cdff14186f1af11fae865a5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>stepsRemaining</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a6d35d4037ec7b5837c8dd84dde6e0264</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt;</type>
      <name>getPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a9c8fa37c580b6e3dc5305ab2269cb25a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt;</type>
      <name>getGriddingAlgorithm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a98f3f1b97d50eada757f0e3f611ba3d7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerVisItDataWriter</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a09f8bc452a3344c9760b9071679f2c95</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::appu::VisItDataWriter&lt; NDIM &gt; &gt; visit_writer)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setupPlotData</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a3df2d7b7bfd2c26cba31e41abca6cd95</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getNumberOfCycles</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ac2849e1c072c230c7a9866078d11c2bf</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getCurrentCycleNumber</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ae67894b4d0915c8b6250706c6b912160</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>getCurrentTimeStepSize</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>aa403af58a2b538a6dc231905c758293a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>preprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ad37decc9a9cb33397c30b1b885324a0c</anchor>
      <arglist>(double current_time, double new_time, int num_cycles=1)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>integrateHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a214c3c38f662cf246f1d20f07c5badc9</anchor>
      <arglist>(double current_time, double new_time, int cycle_num=0)=0</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>skipCycle</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ab05cac4c0946b1864a983e6d1ce26ca3</anchor>
      <arglist>(double current_time, double new_time, int cycle_num=0)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>postprocessIntegrateHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a033261287aebd6ceaaadd8f989c5a525</anchor>
      <arglist>(double current_time, double new_time, bool skip_synchronize_new_state_data, int num_cycles=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerPreprocessIntegrateHierarchyCallback</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ac9229eef598ccd178b5f9c5a49f876cf</anchor>
      <arglist>(PreprocessIntegrateHierarchyCallbackFcnPtr callback, void *ctx=NULL)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerIntegrateHierarchyCallback</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a4ff05eb501095d3f4ee4115c2af1a159</anchor>
      <arglist>(IntegrateHierarchyCallbackFcnPtr callback, void *ctx=NULL)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerPostprocessIntegrateHierarchyCallback</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a0103dd61459dcba961f3fa13d0877742</anchor>
      <arglist>(PostprocessIntegrateHierarchyCallbackFcnPtr callback, void *ctx=NULL)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerApplyGradientDetectorCallback</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a6f8e39591753e05e4242ef27a5bdffa9</anchor>
      <arglist>(ApplyGradientDetectorCallbackFcnPtr callback, void *ctx=NULL)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ae5445debcbe3d436e63b08c20a6572fb</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt;(NULL), bool allocate_data=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetHierarchyConfiguration</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>afc6fbb35718daa8112f5e0fcf4874518</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyGradientDetector</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a29746f4e75b726a9ebc85c533ef8e0cf</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt;</type>
      <name>getContext</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a9bcbbc9f1522182beacba6edb31e747c</anchor>
      <arglist>(VariableContextType ctx_type) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt;</type>
      <name>getCurrentContext</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a5f7635fb69e14e9520d76435111ec527</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt;</type>
      <name>getNewContext</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>affc34b04cbc6c10658b16e3781c16801</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt;</type>
      <name>getScratchContext</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>af01429ba43946b17e304fd4b9151ec88</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isAllocatedPatchData</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>af763a6880584b8dea380c1e63a4f33db</anchor>
      <arglist>(int data_idx, int coarsest_ln=-1, int finest_ln=-1) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>allocatePatchData</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a2c75f2c19af9a09a0ee29554b2288078</anchor>
      <arglist>(int data_idx, double data_time, int coarsest_ln=-1, int finest_ln=-1) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocatePatchData</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a05cec06fd5e990b3b788840cc82c0c22</anchor>
      <arglist>(int data_idx, int coarsest_ln=-1, int finest_ln=-1) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; HierarchyMathOps &gt;</type>
      <name>getHierarchyMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a9316a56b1f5a406ac6348d38e4b25056</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>aa596e4a54c175070d5362607af7b1329</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual double</type>
      <name>getMaximumTimeStepSizeSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ac11a5e92deff799a68aa71dca9027ca3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>synchronizeHierarchyDataSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>acf59436ba46515415117f6cecd1097bc</anchor>
      <arglist>(VariableContextType ctx_type)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>resetTimeDependentHierarchyDataSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a1bfa5de2713a4f7b044e3a7f59950501</anchor>
      <arglist>(double new_time)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>resetIntegratorToPreadvanceStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>af61f6ef496ce51b07947809c46b79c5d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual bool</type>
      <name>atRegridPointSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a7225ceeb7b0a6cecc4bcf428d6e44b08</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>setupPlotDataSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a40ee7f1fa185a4652265735c76ba1113</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeLevelDataSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a2bddf84a64fd161a3b442bd90c01c8cc</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level, bool allocate_data)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>resetHierarchyConfigurationSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a628c9f032715af4f4a506b27d0ac0b7e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyGradientDetectorSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a692855d3982df504232d2dd70a5b46eb</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>putToDatabaseSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a50775b336487c703984b33db64edbb49</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>executePreprocessIntegrateHierarchyCallbackFcns</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a7b75686ce9dd26b40409f734743dbcd1</anchor>
      <arglist>(double current_time, double new_time, int num_cycles)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>executeIntegrateHierarchyCallbackFcns</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>abc2a4a64d6da28445c55d8df0999fa88</anchor>
      <arglist>(double current_time, double new_time, int cycle_num)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>executePostprocessIntegrateHierarchyCallbackFcns</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a02c4cb8da4a402b481a6f8d70101058a</anchor>
      <arglist>(double current_time, double new_time, bool skip_synchronize_new_state_data, int num_cycles)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>executeApplyGradientDetectorCallbackFcns</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a49fd60c3de06c89917322fe2157b0b6e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerVariable</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a2b2ce04d072210438e34aab7de9f97e1</anchor>
      <arglist>(int &amp;current_idx, int &amp;new_idx, int &amp;scratch_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; variable, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;scratch_ghosts=SAMRAI::hier::IntVector&lt; NDIM &gt;(0), const std::string &amp;coarsen_name=&quot;NO_COARSEN&quot;, const std::string &amp;refine_name=&quot;NO_REFINE&quot;, SAMRAI::tbox::Pointer&lt; CartGridFunction &gt; init_fcn=SAMRAI::tbox::Pointer&lt; CartGridFunction &gt;(NULL))</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerVariable</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a2b6db0abbf486c756d708cbbb3e0e335</anchor>
      <arglist>(int &amp;idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; variable, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghosts=SAMRAI::hier::IntVector&lt; NDIM &gt;(0), SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt; ctx=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt;(NULL))</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerGhostfillRefineAlgorithm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a888b46b781fe70cb443583cb41711ee7</anchor>
      <arglist>(const std::string &amp;name, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineAlgorithm&lt; NDIM &gt; &gt; ghostfill_alg, SAMRAI::xfer::RefinePatchStrategy&lt; NDIM &gt; *ghostfill_patch_strategy=NULL)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerProlongRefineAlgorithm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a936f032b0435ca04a7e8fbfaaa6218a8</anchor>
      <arglist>(const std::string &amp;name, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineAlgorithm&lt; NDIM &gt; &gt; prolong_alg, SAMRAI::xfer::RefinePatchStrategy&lt; NDIM &gt; *prolong_patch_strategy=NULL)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerCoarsenAlgorithm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a902fb9dd6a4e4d641f6f2108fd5915fd</anchor>
      <arglist>(const std::string &amp;name, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenAlgorithm&lt; NDIM &gt; &gt; coarsen_alg, SAMRAI::xfer::CoarsenPatchStrategy&lt; NDIM &gt; *coarsen_patch_strategy=NULL)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineAlgorithm&lt; NDIM &gt; &gt;</type>
      <name>getGhostfillRefineAlgorithm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a4fadb95fe1a67aa41d8dcced48241229</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineAlgorithm&lt; NDIM &gt; &gt;</type>
      <name>getProlongRefineAlgorithm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a32a39a6214727ac8c7bba4cdbb0bb23b</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenAlgorithm&lt; NDIM &gt; &gt;</type>
      <name>getCoarsenAlgorithm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a871779a56f7e379bd69c134246c946af</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;</type>
      <name>getGhostfillRefineSchedules</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a9639de1b88f7a364161374b704f5f079</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;</type>
      <name>getProlongRefineSchedules</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>aa2f695110e81eda9afc2a833d80adc7c</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;</type>
      <name>getCoarsenSchedules</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a8857a65cde86f56c046df3aad0408a6c</anchor>
      <arglist>(const std::string &amp;name) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerChildHierarchyIntegrator</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>abd2fd0095daa3b69d308b04408959279</anchor>
      <arglist>(HierarchyIntegrator *child_integrator)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>registerParentHierarchyIntegrator</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>ae62554a0debc2042ce923eaa38997b88</anchor>
      <arglist>(HierarchyIntegrator *parent_integrator)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; HierarchyMathOps &gt;</type>
      <name>buildHierarchyMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a644ed9395afdcf4e7a0a3f2e705e9c01</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>setupTagBuffer</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a9a52774e50c70af5ef1d2f22d6b38420</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; NDIM &gt; &gt; gridding_alg)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>regriddingHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a534118859ede84a764aa3cbc2a4c9532</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>atRegridTimeStep</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a5d2fb76576afb44f7cd822348926eae7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::VariableContext &gt;</type>
      <name>d_current_context</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a4b746c6b3e47b2dd0824d9f980c7a098</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SAMRAI::hier::ComponentSelector</type>
      <name>d_fill_after_regrid_bc_idxs</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a1c39ac5f83a1d5cce38f66927931245a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::vector&lt; PreprocessIntegrateHierarchyCallbackFcnPtr &gt;</type>
      <name>d_preprocess_integrate_hierarchy_callbacks</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a804f74ba0a883d00bacca2ed4dd16c3c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected" static="yes">
      <type>static const std::string</type>
      <name>SYNCH_CURRENT_DATA_ALG</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_integrator.html</anchorfile>
      <anchor>a2a71c9e24e273b7a91e6c27f845c1250</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::HierarchyMathOps</name>
    <filename>class_i_b_t_k_1_1_hierarchy_math_ops.html</filename>
    <member kind="function">
      <type></type>
      <name>HierarchyMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a99d32688cd432b1cd2aed61b009e83c5</anchor>
      <arglist>(const std::string &amp;name, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_ln=-1, int finest_ln=-1, const std::string &amp;coarsen_op_name=&quot;CONSERVATIVE_COARSEN&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~HierarchyMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a3cf0c5a9ee4744a5c4dcf44fbda272b7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a29fc8140fab8b98a3d1689df4b3f0d0c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetLevels</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a3fe79b486c12a8c9b1e292407e425957</anchor>
      <arglist>(int coarsest_ln, int finest_ln)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt;</type>
      <name>getCellWeightVariable</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>ac80fc51b3c48cbf0f80c479082b1fb00</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getCellWeightPatchDescriptorIndex</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>acc263c15026d5879aaafb6ea95a2385f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt;</type>
      <name>getFaceWeightVariable</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>ac2c5f41e47f28bf70c1ab7c09b2b8aec</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getFaceWeightPatchDescriptorIndex</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>ad520f83c3bf095e6731231c4904e5256</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt;</type>
      <name>getSideWeightVariable</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>acb7ec2d48ed15e5d2fbf6926707c7c35</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getSideWeightPatchDescriptorIndex</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>ac3b343a7b7f5ac7a972ee18dde46be4a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getVolumeOfPhysicalDomain</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a631b1ce26adc223b515b337a28ebba4e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCoarsenOperatorName</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a1a5addb7d3f813f07c97e3ed0c6c12f0</anchor>
      <arglist>(const std::string &amp;coarsen_op_name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a7ab7faa43ffa0cdd8b2602cce5ed6572</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a00acb3173655d5571b82aa6daab06040</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a8b4b4903f64f82dcbdeac28a6aa975ba</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a112fb85f11c256ddc9a527327dcae21c</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a5192aa2d731409b0b5a47eca13de6b77</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a6f4770322ef3f1e5d5fb18c40faec758</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a1e766aabef2159c1a9c74b139424a7d2</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rot</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a3e47a236737ef8c952a228ca5c79fd24</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rot</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a70cca384bc4c65514223edb03aa3f946</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rot</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a52fa6d81b715cc5ac5c5cec49e6c853a</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rot</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a7fcf9a5c04493453609393cea2784f61</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>div</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a14fc391101588392a8d7c8f4773c92b8</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, double alpha, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src2_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>div</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>ac93b06312a1828cc54edfd9387f7aaa6</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, double alpha, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, bool src1_cf_bdry_synch, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src2_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>div</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a56343419b83dc1a9f5fdfa22baeb6d9a</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, double alpha, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, bool src1_cf_bdry_synch, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src2_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>acca9290c6b64e871498f44c867430e97</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, double alpha, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int src1_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a37628cd821fef66c4e028a8e393cd0c8</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; dst_var, bool dst_cf_bdry_synch, double alpha, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int src1_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a4fc22298ce1a1113df8b922b6a524291</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, bool dst_cf_bdry_synch, double alpha, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int src1_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>aa02b8010eb39acb02ded472941fbf36f</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int src1_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a1d7cbf2a54e608e7e43b0fe4c4a99464</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int src1_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a692e34ea64f1a75cef318d63997dd6a2</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; dst_var, bool dst_cf_bdry_synch, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int src1_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a175d24d4ab27116abfb0f9d6f4fb43e6</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, bool dst_cf_bdry_synch, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int src1_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a0debb1f3620f3945129dc00bd79c5335</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time, bool src_cf_bdry_synch)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>adba421107047d334ecf8f50823927cdd</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time, bool src_cf_bdry_synch)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a378839bb6b14333500d944e90fca8062</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; dst_var, bool dst_cf_bdry_synch, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>ac14dd1886eff1651f0e8416986aedfbc</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, bool dst_cf_bdry_synch, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src_ghost_fill, double src_ghost_fill_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>laplace</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a1756f7733fe029c170a9591c0b63507d</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, double gamma=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>laplace</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>ac591e2df045aefcb258f352c9c294f41</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, double gamma=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src2_var=NULL)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>vc_laplace</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a96498bc9aebc08af9fbafa24ba1bede1</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, double alpha, double beta, int coef_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; coef_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src1_var, SAMRAI::tbox::Pointer&lt; HierarchyGhostCellInterpolation &gt; src1_ghost_fill, double src1_ghost_fill_time, double gamma=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src2_var=NULL)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a94156efced025822191c7a0c251da26d</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, double alpha, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a1ae11de13f0d449e8a3bec3728aed2ff</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0, int alpha_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>af29f6679b7897a1a2089285a654d7d9e</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src1_var, int beta_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; beta_var=NULL, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0, int alpha_depth=0, int beta_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a01695a3819aa0628885e585b0b4d3da1</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; dst_var, double alpha, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src1_var, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a04fa78ccbbba4c41166d5ec0c2f1ba0d</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; dst_var, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src1_var, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0, int alpha_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>af91a37115446131240ade187e6448c83</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; dst_var, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src1_var, int beta_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; beta_var=NULL, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0, int alpha_depth=0, int beta_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a0efdbe13457f3cf4d83b4f1b347041cd</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; dst_var, double alpha, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; src1_var, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a74e050cb556215efb10f8fd0030cec11</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; dst_var, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; src1_var, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0, int alpha_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a6fcce67493240d0ded4cb76ddd0a82ca</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; dst_var, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; src1_var, int beta_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; beta_var=NULL, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0, int alpha_depth=0, int beta_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>adbd6bf6e2e9c062b22244020c97a6cfb</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, double alpha, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src1_var, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>af4c5d6e0a80e8c68cad1637fd3aad19d</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src1_var, double beta=0.0, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0, int alpha_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>af45edfe5f8533bbc198306318b6c417b</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; dst_var, int alpha_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; alpha_var, int src1_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src1_var, int beta_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; beta_var=NULL, int src2_idx=-1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideVariable&lt; NDIM, double &gt; &gt; src2_var=NULL, int dst_depth=0, int src1_depth=0, int src2_depth=0, int alpha_depth=0, int beta_depth=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseL1Norm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>ab5409cd3d9c6da21be0c7f46854af2b7</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseL2Norm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a080c1a60d98583964bc35f134b29f311</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMaxNorm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a223fa153f15590c457c72ac13284d10b</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellVariable&lt; NDIM, double &gt; &gt; src_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseL1Norm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>aab8d2cd0d3f0f30e8509e2b38650f05a</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; src_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseL2Norm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>a26374e2b2e1843a2b0311985d822989a</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; src_var)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMaxNorm</name>
      <anchorfile>class_i_b_t_k_1_1_hierarchy_math_ops.html</anchorfile>
      <anchor>aa240a28a2274cff2bac4d097481e082a</anchor>
      <arglist>(int dst_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; dst_var, int src_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeVariable&lt; NDIM, double &gt; &gt; src_var)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::IndexUtilities</name>
    <filename>class_i_b_t_k_1_1_index_utilities.html</filename>
    <member kind="function" static="yes">
      <type>static SAMRAI::hier::Index&lt; NDIM &gt;</type>
      <name>getCellIndex</name>
      <anchorfile>class_i_b_t_k_1_1_index_utilities.html</anchorfile>
      <anchor>a14433b3c19ce5073a56ab8b18473172f</anchor>
      <arglist>(const DoubleArray &amp;X, const double *x_lower, const double *x_upper, const double *dx, const SAMRAI::hier::Index&lt; NDIM &gt; &amp;ilower, const SAMRAI::hier::Index&lt; NDIM &gt; &amp;iupper)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::JacobianOperator</name>
    <filename>class_i_b_t_k_1_1_jacobian_operator.html</filename>
    <base>IBTK::LinearOperator</base>
    <member kind="function">
      <type></type>
      <name>JacobianOperator</name>
      <anchorfile>class_i_b_t_k_1_1_jacobian_operator.html</anchorfile>
      <anchor>a7d5cc63072a95ab18d01c5b9a8688d71</anchor>
      <arglist>(const std::string &amp;object_name)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~JacobianOperator</name>
      <anchorfile>class_i_b_t_k_1_1_jacobian_operator.html</anchorfile>
      <anchor>a571e90b327397703d140abeed9d9eca3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>formJacobian</name>
      <anchorfile>class_i_b_t_k_1_1_jacobian_operator.html</anchorfile>
      <anchor>a5cdc5dd2dcc69a063856e0dcf9e16018</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt;</type>
      <name>getBaseVector</name>
      <anchorfile>class_i_b_t_k_1_1_jacobian_operator.html</anchorfile>
      <anchor>abd460f3b2f27080dd46330b2fe9a7dee</anchor>
      <arglist>() const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::KrylovLinearSolver</name>
    <filename>class_i_b_t_k_1_1_krylov_linear_solver.html</filename>
    <base>IBTK::LinearSolver</base>
    <member kind="function">
      <type></type>
      <name>KrylovLinearSolver</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver.html</anchorfile>
      <anchor>ab85c7bafc1fe8a48e14b23381f28783b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~KrylovLinearSolver</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver.html</anchorfile>
      <anchor>af6b5339bbd19582d575c06a123aa02c0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHierarchyMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver.html</anchorfile>
      <anchor>aad431cae5aaa8b41c60868df94186fe9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; HierarchyMathOps &gt; hier_math_ops)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver.html</anchorfile>
      <anchor>a5f15a8d19ef4e7eeb34a872156af9e06</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSolutionTime</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver.html</anchorfile>
      <anchor>abcfb2eed8c264ab0d7be76698de72eeb</anchor>
      <arglist>(double solution_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver.html</anchorfile>
      <anchor>a081bc2d246d5303c021f33cfbe309271</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setOperator</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver.html</anchorfile>
      <anchor>a54f3c32afdd1a452fd0e4a072411efe8</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LinearOperator &gt; A)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual SAMRAI::tbox::Pointer&lt; LinearOperator &gt;</type>
      <name>getOperator</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver.html</anchorfile>
      <anchor>a9d79dca6a53fa3cd791f98c24968340e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver.html</anchorfile>
      <anchor>a427ec5564a20480e0c3a79311585843b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LinearSolver &gt; pc_solver=NULL)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual SAMRAI::tbox::Pointer&lt; LinearSolver &gt;</type>
      <name>getPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver.html</anchorfile>
      <anchor>a4f058f0217fba01f9718709aae76c76b</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::KrylovLinearSolverManager</name>
    <filename>class_i_b_t_k_1_1_krylov_linear_solver_manager.html</filename>
    <member kind="typedef">
      <type>SAMRAI::tbox::Pointer&lt; KrylovLinearSolver &gt;(*</type>
      <name>SolverMaker</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_manager.html</anchorfile>
      <anchor>af03c7b4fe6e84d29abce1b33dba4256c</anchor>
      <arglist>)(const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; KrylovLinearSolver &gt;</type>
      <name>allocateSolver</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_manager.html</anchorfile>
      <anchor>a195878c1498c0725146bcf16b3e97ed8</anchor>
      <arglist>(const std::string &amp;solver_type, const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerSolverFactoryFunction</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_manager.html</anchorfile>
      <anchor>a727ceb8b9408613ebd69472a64ae5bd5</anchor>
      <arglist>(const std::string &amp;solver_type, SolverMaker solver_maker)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static KrylovLinearSolverManager *</type>
      <name>getManager</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_manager.html</anchorfile>
      <anchor>a203f6c91ca1cdec764743aa4f1f8ff07</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>freeManager</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_manager.html</anchorfile>
      <anchor>af907f1f366e758deae39b813ebfebc7f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>UNDEFINED</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_manager.html</anchorfile>
      <anchor>a1d782b9f98a30d3815e047709f2ef57a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_manager.html</anchorfile>
      <anchor>a277097208c9740fa27703597c438cb5d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>KrylovLinearSolverManager</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_manager.html</anchorfile>
      <anchor>a771b38db50af2d1490fc4a81db01d34a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>~KrylovLinearSolverManager</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_manager.html</anchorfile>
      <anchor>a5e8411665e91119d66215b3499e0f61b</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::KrylovLinearSolverPoissonSolverInterface</name>
    <filename>class_i_b_t_k_1_1_krylov_linear_solver_poisson_solver_interface.html</filename>
    <base>IBTK::PoissonSolver</base>
    <member kind="function">
      <type></type>
      <name>KrylovLinearSolverPoissonSolverInterface</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_poisson_solver_interface.html</anchorfile>
      <anchor>a65139838e247dda6c660bcaf908d122c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~KrylovLinearSolverPoissonSolverInterface</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_poisson_solver_interface.html</anchorfile>
      <anchor>a455f9b44c4a7977700c4c5e0f5139221</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPoissonSpecifications</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_poisson_solver_interface.html</anchorfile>
      <anchor>ada19e41118d66f61636a8abd100b29e9</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoef</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_poisson_solver_interface.html</anchorfile>
      <anchor>a636af79c413694a735d4ec18fdd0d6b6</anchor>
      <arglist>(SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_t_k_1_1_krylov_linear_solver_poisson_solver_interface.html</anchorfile>
      <anchor>ab808a58d9d0d42d749638e322a664085</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LaplaceOperator</name>
    <filename>class_i_b_t_k_1_1_laplace_operator.html</filename>
    <base>IBTK::LinearOperator</base>
  </compound>
  <compound kind="class">
    <name>IBTK::LData</name>
    <filename>class_i_b_t_k_1_1_l_data.html</filename>
    <base>SAMRAI::tbox::Serializable</base>
    <member kind="function">
      <type></type>
      <name>LData</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a308cb60cc812943a09c32c36e263971f</anchor>
      <arglist>(const std::string &amp;name, unsigned int num_local_nodes, unsigned int depth, const std::vector&lt; int &gt; &amp;nonlocal_petsc_indices=std::vector&lt; int &gt;(0))</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LData</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a18401ad8f94f944fbec772fe1ae4f285</anchor>
      <arglist>(const std::string &amp;name, Vec vec, const std::vector&lt; int &gt; &amp;nonlocal_petsc_indices=std::vector&lt; int &gt;(0), const bool manage_petsc_vec=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LData</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a1f6f158b5677963e1e9338a4ee6be4bd</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LData</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a1462aea8515aa41442a40722ed9ac078</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>resetData</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a8ad4fea9a202cacb918eeeca4bdd97ec</anchor>
      <arglist>(Vec vec, const std::vector&lt; int &gt; &amp;nonlocal_petsc_indices=std::vector&lt; int &gt;(0), const bool manage_petsc_vec=true)</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getName</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a6eff55128ee78646d1e044783071dfa9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getGlobalNodeCount</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a745e79c5ee9205618aa1a5e1644ec6a1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getLocalNodeCount</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>ab70a6ee73a97af86a62ed60190cbf521</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getGhostNodeCount</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a2ec35bdca9a59464eb1b0823ecb2b28e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getDepth</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a73f733c1ddaeb20b1bc48ae97e812844</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Vec</type>
      <name>getVec</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>ae837fb356ffaf419beac89314348bcb3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>boost::multi_array_ref&lt; double, 1 &gt; *</type>
      <name>getArray</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a162aaaa21135a8fd850207db97783599</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>boost::multi_array_ref&lt; double, 1 &gt; *</type>
      <name>getLocalFormArray</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a0d74b27ea9d3b49a8dfba3e843f6f5af</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>boost::multi_array_ref&lt; double, 1 &gt; *</type>
      <name>getGhostedLocalFormArray</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>ac06409d2b3cb7490d0aeeeec991821f5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>boost::multi_array_ref&lt; double, 2 &gt; *</type>
      <name>getVecArray</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a55cec55c85a3bd877f54bfb14b9ef6cc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>boost::multi_array_ref&lt; double, 2 &gt; *</type>
      <name>getLocalFormVecArray</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>aeecc2236a25bbb816b35b0650c6fe96d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>boost::multi_array_ref&lt; double, 2 &gt; *</type>
      <name>getGhostedLocalFormVecArray</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a41de3424869163ff8fe311b209893721</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>restoreArrays</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>aa78d834678ec57baea3fc552bc2bc3c0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>beginGhostUpdate</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>aa6fcf7727e8048046c555d1cb694821f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>endGhostUpdate</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>ada12ad1d2d027884e591d1657c1a5d3e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_t_k_1_1_l_data.html</anchorfile>
      <anchor>a2a7afc3fa8f430e98798e977056e7757</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LDataManager</name>
    <filename>class_i_b_t_k_1_1_l_data_manager.html</filename>
    <base>SAMRAI::tbox::Serializable</base>
    <base>StandardTagAndInitStrategy&lt; NDIM &gt;</base>
    <member kind="function">
      <type>const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getGhostCellWidth</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a15ac4f21e711c807aaa4de9ed4ec9b51</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getDefaultInterpKernelFunction</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ac4959bb0c248316fd90fcfc6cbaa3f3d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getDefaultSpreadKernelFunction</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>afef6e0e6c127a8e2031728713893c26f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a84a7a199cfd186f610f87c5117e699df</anchor>
      <arglist>(int f_data_idx, SAMRAI::tbox::Pointer&lt; LData &gt; F_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, SAMRAI::tbox::Pointer&lt; LData &gt; ds_data, RobinPhysBdryPatchStrategy *f_phys_bdry_op, int level_num, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0, bool F_data_ghost_node_update=true, bool X_data_ghost_node_update=true, bool ds_data_ghost_node_update=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a8bec837fe1c813bd4b51829f574de630</anchor>
      <arglist>(int f_data_idx, SAMRAI::tbox::Pointer&lt; LData &gt; F_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, SAMRAI::tbox::Pointer&lt; LData &gt; ds_data, const std::string &amp;spread_kernel_fcn, RobinPhysBdryPatchStrategy *f_phys_bdry_op, int level_num, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0, bool F_data_ghost_node_update=true, bool X_data_ghost_node_update=true, bool ds_data_ghost_node_update=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a978d332cf68b579d81dd176708798d29</anchor>
      <arglist>(int f_data_idx, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;F_data, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;X_data, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;ds_data, RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0, bool F_data_ghost_node_update=true, bool X_data_ghost_node_update=true, bool ds_data_ghost_node_update=true, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>aa1c2d501c1459852b9d57faa00f9a1c2</anchor>
      <arglist>(int f_data_idx, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;F_data, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;X_data, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;ds_data, const std::string &amp;spread_kernel_fcn, RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0, bool F_data_ghost_node_update=true, bool X_data_ghost_node_update=true, bool ds_data_ghost_node_update=true, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>abfe90013688c6dbaf9d827cc652234f0</anchor>
      <arglist>(int f_data_idx, SAMRAI::tbox::Pointer&lt; LData &gt; F_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, RobinPhysBdryPatchStrategy *f_phys_bdry_op, int level_num, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0, bool F_data_ghost_node_update=true, bool X_data_ghost_node_update=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>af3b321e1f75f119445da3a107d43cd76</anchor>
      <arglist>(int f_data_idx, SAMRAI::tbox::Pointer&lt; LData &gt; F_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, const std::string &amp;spread_kernel_fcn, RobinPhysBdryPatchStrategy *f_phys_bdry_op, int level_num, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0, bool F_data_ghost_node_update=true, bool X_data_ghost_node_update=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a61fa9b0dda27aedcf29435addd0d1ad1</anchor>
      <arglist>(int f_data_idx, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;F_data, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;X_data, RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0, bool F_data_ghost_node_update=true, bool X_data_ghost_node_update=true, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a8f543f7890526c6b6f2e4d6466cce65e</anchor>
      <arglist>(int f_data_idx, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;F_data, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;X_data, const std::string &amp;spread_kernel_fcn, RobinPhysBdryPatchStrategy *f_phys_bdry_op, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_prolongation_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0, bool F_data_ghost_node_update=true, bool X_data_ghost_node_update=true, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ab4648251121262ca5bae15054fb760a6</anchor>
      <arglist>(int f_data_idx, SAMRAI::tbox::Pointer&lt; LData &gt; F_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, int level_num, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_synch_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt;(), const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_ghost_fill_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ae93643716df661265a436ba6660e9601</anchor>
      <arglist>(int f_data_idx, SAMRAI::tbox::Pointer&lt; LData &gt; F_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, const std::string &amp;interp_kernel_fcn, int level_num, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_synch_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt;(), const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_ghost_fill_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a9d16ba694da412559a82aa55720c76ce</anchor>
      <arglist>(int f_data_idx, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;F_data, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;X_data, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_synch_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt;(), const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_ghost_fill_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ae492e87313171335aa52e59050c88091</anchor>
      <arglist>(int f_data_idx, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;F_data, std::vector&lt; SAMRAI::tbox::Pointer&lt; LData &gt; &gt; &amp;X_data, const std::string &amp;interp_kernel_fcn, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_synch_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::CoarsenSchedule&lt; NDIM &gt; &gt; &gt;(), const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt; &amp;f_ghost_fill_scheds=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; &gt;(), double fill_data_time=0.0, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLInitStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a56b81e163585cde632beaea203a574a1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LInitStrategy &gt; lag_init)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>freeLInitStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a0c54478bd27a46546e19292ed3f2ef08</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerVisItDataWriter</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>aa68d470c1bb587fb77525d7a303bedcd</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::appu::VisItDataWriter&lt; NDIM &gt; &gt; visit_writer)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLSiloDataWriter</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>aa9946bcc90dc2bd5841fca70b8fa82a2</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LSiloDataWriter &gt; silo_writer)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLoadBalancer</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a4091bdb848f92c5ee87d5c6eb92145d0</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::mesh::LoadBalancer&lt; NDIM &gt; &gt; load_balancer, int workload_data_idx)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>levelContainsLagrangianData</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ad0c98065e53101d60eb41cc0a432ab36</anchor>
      <arglist>(int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getNumberOfNodes</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a51dab90c236af236ce6d67065e3a8dd5</anchor>
      <arglist>(int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getNumberOfLocalNodes</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a0b2fae9e1759351c9ad1c10cf8f6a5bf</anchor>
      <arglist>(int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getNumberOfGhostNodes</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a9469a6622ec345fc726cf9cb04dcaca1</anchor>
      <arglist>(int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getGlobalNodeOffset</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a4c4daacc4e65eb79dcc4149e05e4180e</anchor>
      <arglist>(int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; LMesh &gt;</type>
      <name>getLMesh</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a907cfd523fc51dff6915639e1e70dcef</anchor>
      <arglist>(int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; LData &gt;</type>
      <name>getLData</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a9e527a4f0d56cd7f8197cbd751c565c7</anchor>
      <arglist>(const std::string &amp;quantity_name, int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; LData &gt;</type>
      <name>createLData</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ae201fcac2bb4b6b3012215e1fbb2b0fa</anchor>
      <arglist>(const std::string &amp;quantity_name, int level_number, unsigned int depth=1, bool maintain_data=false)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getLNodePatchDescriptorIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>aa5c728b2718eb824944b0a903f034bc0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getWorkloadPatchDescriptorIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ac2025d2865d330a55cd86a204d44cfc1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNodeCountPatchDescriptorIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ad8c576c8f9d4e3ddea088ab6012bf561</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::string &gt;</type>
      <name>getLagrangianStructureNames</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a7db82e1e7d8ce9787d0559769f760c0a</anchor>
      <arglist>(int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; int &gt;</type>
      <name>getLagrangianStructureIDs</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a44433e232b65cc0594258f93c340a23f</anchor>
      <arglist>(int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getLagrangianStructureID</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a05422c8e5b41c764ac86b2aae9e27d3c</anchor>
      <arglist>(int lagrangian_index, int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getLagrangianStructureID</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a31e75c1a713497463ea2ea32ea12a6c5</anchor>
      <arglist>(const std::string &amp;structure_name, int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getLagrangianStructureName</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a66d45be24fe35aeb3bf4d9d75ec5edee</anchor>
      <arglist>(int structure_id, int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; int, int &gt;</type>
      <name>getLagrangianStructureIndexRange</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ada6587e55e27de16ea87b8b1a67fcdeb</anchor>
      <arglist>(int structure_id, int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>Point</type>
      <name>computeLagrangianStructureCenterOfMass</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a6d2d121ec746f96f49662e4e6478b86a</anchor>
      <arglist>(int structure_id, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; Point, Point &gt;</type>
      <name>computeLagrangianStructureBoundingBox</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>abe056cccd40e450593e73f0342a3c297</anchor>
      <arglist>(int structure_id, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinitLagrangianStructure</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>afdaafde7e8e58a5f2ab5ee5deaeb5ca0</anchor>
      <arglist>(const Point &amp;X_center, int structure_id, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>displaceLagrangianStructure</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a11e91185ad0ce241dafd8fa176dea8fc</anchor>
      <arglist>(const Vector &amp;dX, int structure_id, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>activateLagrangianStructures</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>aef4eaedcdd89b31a5dfb6a668fbb3807</anchor>
      <arglist>(const std::vector&lt; int &gt; &amp;structure_ids, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>inactivateLagrangianStructures</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a3aeb7b0d36f493ac4e242da1f73b5598</anchor>
      <arglist>(const std::vector&lt; int &gt; &amp;structure_ids, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>getLagrangianStructureIsActivated</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a8891767ebaaab58b10e0667afd974788</anchor>
      <arglist>(int structure_id, int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>zeroInactivatedComponents</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a2a53da32846c4e28e4caf63e5e939bd7</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LData &gt; lag_data, int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>mapLagrangianToPETSc</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a1743485b455af4a105bf193a93478e52</anchor>
      <arglist>(std::vector&lt; int &gt; &amp;inds, int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>mapPETScToLagrangian</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ac4f25818730b643481c00c5c3bf56417</anchor>
      <arglist>(std::vector&lt; int &gt; &amp;inds, int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>scatterLagrangianToPETSc</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a2dfc2cbadd6c8d4c705640e77e5e5666</anchor>
      <arglist>(Vec &amp;lagrangian_vec, Vec &amp;petsc_vec, int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>scatterPETScToLagrangian</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a31744cbe297ffea1cbd1cb80d4ab712e</anchor>
      <arglist>(Vec &amp;petsc_vec, Vec &amp;lagrangian_vec, int level_number) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>scatterToAll</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a44eeb7c130263f1a8d363a8656231887</anchor>
      <arglist>(Vec &amp;parallel_vec, Vec &amp;sequential_vec) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>scatterToZero</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ae93c3331c8523a28590196ab74b77810</anchor>
      <arglist>(Vec &amp;parallel_vec, Vec &amp;sequential_vec) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>beginDataRedistribution</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a928620cb87d696a8584728f6d4ff4568</anchor>
      <arglist>(int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>endDataRedistribution</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a03ef7c59981d1ba33d9fb14cc13d741c</anchor>
      <arglist>(int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>updateWorkloadEstimates</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a56c9c960bf0e7ce517474eef4849cf7b</anchor>
      <arglist>(int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>updateNodeCountData</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a43c25f4a3ce1c1e872ce214766db1aba</anchor>
      <arglist>(int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>acaee193e580c056f4473af21590ef333</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt;(NULL), bool allocate_data=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetHierarchyConfiguration</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a05ec58dc8cedf01abb8b5b9a56ff776a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_ln, int finest_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyGradientDetector</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>acc21584585cefe0192578e4e68472b87</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a12305cb4e98703b734af75f5475bdcb4</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a36cecd3accdf64b12e8f1093fb3ac7d8</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt;</type>
      <name>getPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a46e4035f89fffe6d3de3d90f31883fe9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchLevels</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>ab66d2c666c08047508864acd1d12f27c</anchor>
      <arglist>(int coarsest_ln, int finest_ln)</arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; int, int &gt;</type>
      <name>getPatchLevels</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>aa2ca4ced55f6fc4b0740dce2bef214b3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" static="yes">
      <type>static LDataManager *</type>
      <name>getManager</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a5e4bfd3d20e12bcd1d5324fc01babbe6</anchor>
      <arglist>(const std::string &amp;name, const std::string &amp;default_interp_kernel_fcn, const std::string &amp;default_spread_kernel_fcn, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;min_ghost_width=SAMRAI::hier::IntVector&lt; NDIM &gt;(0), bool register_for_restart=true)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>freeAllManagers</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a9c93e1e915c0bef2764cfbb62f51e592</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>POSN_DATA_NAME</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a359fffa7d2f68d6ec3ad02d7b8599b0d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>INIT_POSN_DATA_NAME</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>afd892b71bce2c9096c49efd80db288de</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>VEL_DATA_NAME</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a04e7f66b3437a4b69237fef346348ad6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>LDataManager</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>afdaa0f15640377c8f2c46a1d36ad4059</anchor>
      <arglist>(const std::string &amp;object_name, const std::string &amp;default_interp_kernel_fcn, const std::string &amp;default_spread_kernel_fcn, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>~LDataManager</name>
      <anchorfile>class_i_b_t_k_1_1_l_data_manager.html</anchorfile>
      <anchor>a43b47ef5a2be4290da1f1ecd31d956e2</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LEInteractor</name>
    <filename>class_i_b_t_k_1_1_l_e_interactor.html</filename>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>setFromDatabase</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a70a80c371110e56e367a6def06323700</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>printClassData</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>af5ec6b7bb60f6d91b1e009140aa78293</anchor>
      <arglist>(std::ostream &amp;os)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>getStencilSize</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>abc124642b2ff9d49690305c14baa52f6</anchor>
      <arglist>(const std::string &amp;kernel_fcn)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>getMinimumGhostWidth</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a1a497a047fdc0c428d6dc8f92038d74e</anchor>
      <arglist>(const std::string &amp;kernel_fcn)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>ab50bf16bbf5fbb976e1309520104bf4c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LData &gt; Q_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>ab69780eb2132908f84c93bf486697593</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LData &gt; Q_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a07917f8534db959946600bbb089a9535</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LData &gt; Q_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a5fb7d67b996b96b010ae4a1b02aa290b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LData &gt; Q_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a6b49c4e0e351873b1c43077ecbfd1ba3</anchor>
      <arglist>(double *Q_data, int Q_depth, const double *X_data, int X_depth, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>ad9c94e3a75d80a301f3745f6764e434d</anchor>
      <arglist>(double *Q_data, int Q_depth, const double *X_data, int X_depth, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>aeed991dacd0505074bbea8df88a6a1f6</anchor>
      <arglist>(double *Q_data, int Q_depth, const double *X_data, int X_depth, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a3d5b55d76e65b54d27b0cab3b6fb34a7</anchor>
      <arglist>(double *Q_data, int Q_depth, const double *X_data, int X_depth, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a51da1be138dbde8656fd1d9e6221ff7e</anchor>
      <arglist>(std::vector&lt; double &gt; &amp;Q_data, int Q_depth, const std::vector&lt; double &gt; &amp;X_data, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a8fb4b8f4b0d05695b8a25cfce5c6254d</anchor>
      <arglist>(std::vector&lt; double &gt; &amp;Q_data, int Q_depth, const std::vector&lt; double &gt; &amp;X_data, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a471a9f74717fdbb5405cf5c8e3ad584b</anchor>
      <arglist>(std::vector&lt; double &gt; &amp;Q_data, int Q_depth, const std::vector&lt; double &gt; &amp;X_data, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>ad87fd34d2ccac7f0bdf1a1aa49bd7bbd</anchor>
      <arglist>(std::vector&lt; double &gt; &amp;Q_data, int Q_depth, const std::vector&lt; double &gt; &amp;X_data, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a41af71dd97bb361c95933c423e39fda3</anchor>
      <arglist>(double *Q_data, int Q_size, int Q_depth, const double *X_data, int X_size, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>ae28770a31eb0fdfcfe6e29b27fc15247</anchor>
      <arglist>(double *Q_data, int Q_size, int Q_depth, const double *X_data, int X_size, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>ada6962fb97559fd9e710da7bec142332</anchor>
      <arglist>(double *Q_data, int Q_size, int Q_depth, const double *X_data, int X_size, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>interpolate</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a6acc2a62a8dba2d7a7c00ae926c869ed</anchor>
      <arglist>(double *Q_data, int Q_size, int Q_depth, const double *X_data, int X_size, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;interp_box, const std::string &amp;interp_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>aac62e596f20002848dbd03a0113d67f3</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; LData &gt; Q_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a77cc29c3c81d5bb971ee7c07e9feb0bb</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; LData &gt; Q_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a6fc2a84f5731da4c499371444ae6cc79</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; LData &gt; Q_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>af77157126673d3d823fcaa74715931db</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeData&lt; NDIM, double &gt; &gt; q_data, SAMRAI::tbox::Pointer&lt; LData &gt; Q_data, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a1c54a2ec3553c97262a649d56c5eb518</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; q_data, const double *Q_data, int Q_depth, const double *X_data, int X_depth, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a3847e0b185eb578a865525b66fca02e7</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; q_data, const double *Q_data, int Q_depth, const double *X_data, int X_depth, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>aa856f0335e7a9299a3c729d14815d473</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; q_data, const double *Q_data, int Q_depth, const double *X_data, int X_depth, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a28a92589f4f523d397364df5d2ed484c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeData&lt; NDIM, double &gt; &gt; q_data, const double *Q_data, int Q_depth, const double *X_data, int X_depth, SAMRAI::tbox::Pointer&lt; LIndexSetData&lt; T &gt; &gt; idx_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a40563d9c15630b4f760eb7c656d8dbe1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; q_data, const std::vector&lt; double &gt; &amp;Q_data, int Q_depth, const std::vector&lt; double &gt; &amp;X_data, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>ad82cb195990d3d50b84df0fbce6837cb</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; q_data, const std::vector&lt; double &gt; &amp;Q_data, int Q_depth, const std::vector&lt; double &gt; &amp;X_data, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a4ad379896ccebc4acd59428b40d6316e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; q_data, const std::vector&lt; double &gt; &amp;Q_data, int Q_depth, const std::vector&lt; double &gt; &amp;X_data, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a5f26eb425f7e15f7289226992912abb2</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeData&lt; NDIM, double &gt; &gt; q_data, const std::vector&lt; double &gt; &amp;Q_data, int Q_depth, const std::vector&lt; double &gt; &amp;X_data, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a130b9c5a72cfbe9c3dae195519c05f8f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; q_data, const double *Q_data, int Q_size, int Q_depth, const double *X_data, int X_size, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>ac3e53b0cc5869e90e97431167f9277ce</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; q_data, const double *Q_data, int Q_size, int Q_depth, const double *X_data, int X_size, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>ad6e126e69aacb0fad306a1460d895fb4</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; q_data, const double *Q_data, int Q_size, int Q_depth, const double *X_data, int X_size, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>spread</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a1d957ef2085c8b87089fe7c661097ae5</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeData&lt; NDIM, double &gt; &gt; q_data, const double *Q_data, int Q_size, int Q_depth, const double *X_data, int X_size, int X_depth, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;spread_box, const std::string &amp;spread_fcn=&quot;IB_4&quot;)</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static double(*</type>
      <name>s_kernel_fcn</name>
      <anchorfile>class_i_b_t_k_1_1_l_e_interactor.html</anchorfile>
      <anchor>a5b5d8250c83d48ad1fba461f46de5d6b</anchor>
      <arglist>)(double r)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LIndexSetData</name>
    <filename>class_i_b_t_k_1_1_l_index_set_data.html</filename>
    <templarg></templarg>
    <base>IBTK::LSetData</base>
    <member kind="function">
      <type></type>
      <name>LIndexSetData</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>a1baaa6da7ac9cee383e0990c59b6836e</anchor>
      <arglist>(const SAMRAI::hier::Box&lt; NDIM &gt; &amp;box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghosts)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LIndexSetData</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>a8cd7df01a3aa9aa3caffab7c6e06b53b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cacheLocalIndices</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>aa8ba882fc2f6176ad92d900e2d9f9259</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_shift)</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getLagrangianIndices</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>af03194f52d8cb8078ad98fd4c80ae17d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getInteriorLagrangianIndices</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>acfde335fd3f3fa5b8c970152cc059689</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getGhostLagrangianIndices</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>a370c581ec72fac2bbe9d30bf38dcdd14</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getGlobalPETScIndices</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>a1aca997adb9bf245ccdc0971cadbc63f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getInteriorGlobalPETScIndices</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>a8613efc2fbbff01a8c5f1c7d15ec6df9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getGhostGlobalPETScIndices</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>ae35957d86c42ed267dcd32a10f0ae8ff</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getLocalPETScIndices</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>a665500c96d2745a8174d66dec375fbde</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getInteriorLocalPETScIndices</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>a2ea612534167f63d60c01deaffbf6c86</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; int &gt; &amp;</type>
      <name>getGhostLocalPETScIndices</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>a4ca06879a5d432753e2fc7798f4411ae</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; double &gt; &amp;</type>
      <name>getPeriodicShifts</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>afafe67fe390ddf676dd160e97152da10</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; double &gt; &amp;</type>
      <name>getInteriorPeriodicShifts</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>af0199b0f73fe66883bb8674cea080c2e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; double &gt; &amp;</type>
      <name>getGhostPeriodicShifts</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data.html</anchorfile>
      <anchor>a5905753d8a580da800426d11a07541e4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DataIterator</type>
      <name>data_begin</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data.html</anchorfile>
      <anchor>a96c101b247d049f142c57090f0c4f2b6</anchor>
      <arglist>(const SAMRAI::hier::Box&lt; NDIM &gt; &amp;box)</arglist>
    </member>
    <member kind="function">
      <type>DataIterator</type>
      <name>data_end</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data.html</anchorfile>
      <anchor>af6aa59197cdcdd4334af7b56ec9a7491</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LSetData</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data.html</anchorfile>
      <anchor>aeea928dfc67c9ee52209a0fc68081a31</anchor>
      <arglist>(const SAMRAI::hier::Box&lt; NDIM &gt; &amp;box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghosts)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LSetData</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data.html</anchorfile>
      <anchor>aaa73da8b32cc239948d32472e088d12b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="typedef">
      <type>SAMRAI::pdat::CellIterator&lt; NDIM &gt;</type>
      <name>CellIterator</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data.html</anchorfile>
      <anchor>a222852045813115545a48e730fabca4e</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>SAMRAI::pdat::IndexIterator&lt; NDIM, LSet&lt; T &gt;, SAMRAI::pdat::CellGeometry&lt; NDIM &gt; &gt;</type>
      <name>SetIterator</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data.html</anchorfile>
      <anchor>a110476abceb244615b1a0e7085c7fdd2</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>IBTK::LSetDataIterator&lt; T &gt;</type>
      <name>DataIterator</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data.html</anchorfile>
      <anchor>a50ee8be6938b4ff6a8d0db1934d91604</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LIndexSetDataFactory</name>
    <filename>class_i_b_t_k_1_1_l_index_set_data_factory.html</filename>
    <templarg></templarg>
    <base>IBTK::LSetDataFactory</base>
    <member kind="function">
      <type></type>
      <name>LIndexSetDataFactory</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data_factory.html</anchorfile>
      <anchor>affea9cc99fbcdffbd638a8612f458369</anchor>
      <arglist>(const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghosts)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LIndexSetDataFactory</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data_factory.html</anchorfile>
      <anchor>afaac506f7fa405ecd65e7de5369471c3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchData&lt; NDIM &gt; &gt;</type>
      <name>allocate</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data_factory.html</anchorfile>
      <anchor>a382b0e7534aedff04b9f1d806378ada4</anchor>
      <arglist>(const SAMRAI::hier::Box&lt; NDIM &gt; &amp;box, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Arena &gt; pool=NULL) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchData&lt; NDIM &gt; &gt;</type>
      <name>allocate</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data_factory.html</anchorfile>
      <anchor>aa5799c03430e07aa59c8b9c02975187b</anchor>
      <arglist>(const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Arena &gt; pool=NULL) const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSizeOfMemory</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data_factory.html</anchorfile>
      <anchor>a7f3ae48dbef40bc87a7c054bf0d3a9bb</anchor>
      <arglist>(const SAMRAI::hier::Box&lt; NDIM &gt; &amp;box) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchDataFactory&lt; NDIM &gt; &gt;</type>
      <name>cloneFactory</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data_factory.html</anchorfile>
      <anchor>aa969819901a083b8003c5c9eb82c8d13</anchor>
      <arglist>(const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghosts)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>validCopyTo</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_data_factory.html</anchorfile>
      <anchor>afb9caa8ef05c71801173924e0bb942ac</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchDataFactory&lt; NDIM &gt; &gt; &amp;dst_pdf) const </arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LSetDataFactory</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_factory.html</anchorfile>
      <anchor>adf431a57fb70fbaf53ca6d11b43abb78</anchor>
      <arglist>(const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghosts)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LSetDataFactory</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_factory.html</anchorfile>
      <anchor>a6a00a89dcbdf40f1b4ae0c00a3e64e5c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchData&lt; NDIM &gt; &gt;</type>
      <name>allocate</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_factory.html</anchorfile>
      <anchor>a209daf51155574723277d258df4fa4f7</anchor>
      <arglist>(const SAMRAI::hier::Box&lt; NDIM &gt; &amp;box, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Arena &gt; pool=NULL) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchData&lt; NDIM &gt; &gt;</type>
      <name>allocate</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_factory.html</anchorfile>
      <anchor>a3dbcc9866fd167ef5f591601bd137363</anchor>
      <arglist>(const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Arena &gt; pool=NULL) const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSizeOfMemory</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_factory.html</anchorfile>
      <anchor>a98fab4a59af53f7130222024fbc054c8</anchor>
      <arglist>(const SAMRAI::hier::Box&lt; NDIM &gt; &amp;box) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchDataFactory&lt; NDIM &gt; &gt;</type>
      <name>cloneFactory</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_factory.html</anchorfile>
      <anchor>a914a2fdd78c9228b676d86b4c0e0f4ae</anchor>
      <arglist>(const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghosts)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>validCopyTo</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_factory.html</anchorfile>
      <anchor>a3c3aed8682af68d334239033f772171f</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchDataFactory&lt; NDIM &gt; &gt; &amp;dst_pdf) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LIndexSetVariable</name>
    <filename>class_i_b_t_k_1_1_l_index_set_variable.html</filename>
    <templarg></templarg>
    <base>Variable&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>LIndexSetVariable</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_variable.html</anchorfile>
      <anchor>a16334a436db72f50ea580ade41547958</anchor>
      <arglist>(const std::string &amp;name)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LIndexSetVariable</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_variable.html</anchorfile>
      <anchor>a42c0eacc023ab343cbf7bcfc23274e82</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dataLivesOnPatchBorder</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_variable.html</anchorfile>
      <anchor>a27c2102f467e06bd5378801dc6a481c1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>fineBoundaryRepresentsVariable</name>
      <anchorfile>class_i_b_t_k_1_1_l_index_set_variable.html</anchorfile>
      <anchor>a333a01cc3014b5bcce2654e226b0e9ab</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LinearOperator</name>
    <filename>class_i_b_t_k_1_1_linear_operator.html</filename>
    <base>IBTK::GeneralOperator</base>
  </compound>
  <compound kind="class">
    <name>IBTK::LinearSolver</name>
    <filename>class_i_b_t_k_1_1_linear_solver.html</filename>
    <base virtualness="virtual">IBTK::GeneralSolver</base>
  </compound>
  <compound kind="class">
    <name>IBTK::LInitStrategy</name>
    <filename>class_i_b_t_k_1_1_l_init_strategy.html</filename>
    <member kind="function">
      <type></type>
      <name>LInitStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_l_init_strategy.html</anchorfile>
      <anchor>af104c17ab3a990f9c1e2ced49c34377e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LInitStrategy</name>
      <anchorfile>class_i_b_t_k_1_1_l_init_strategy.html</anchorfile>
      <anchor>ae8ac2384f83c9540b0aa8aebc4f69648</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>getLevelHasLagrangianData</name>
      <anchorfile>class_i_b_t_k_1_1_l_init_strategy.html</anchorfile>
      <anchor>a0d72aaee018bb5fe63f08d566341b7c0</anchor>
      <arglist>(int level_number, bool can_be_refined) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual unsigned int</type>
      <name>computeLocalNodeCountOnPatchLevel</name>
      <anchorfile>class_i_b_t_k_1_1_l_init_strategy.html</anchorfile>
      <anchor>a3fb4033ef7503b8f48d671f6314b8277</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeStructureIndexingOnPatchLevel</name>
      <anchorfile>class_i_b_t_k_1_1_l_init_strategy.html</anchorfile>
      <anchor>a257a8ba4f293a6337f9f6ac59bb26089</anchor>
      <arglist>(std::map&lt; int, std::string &gt; &amp;strct_id_to_strct_name_map, std::map&lt; int, std::pair&lt; int, int &gt; &gt; &amp;strct_id_to_lag_idx_range_map, int level_number, double init_data_time, bool can_be_refined, bool initial_time, LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual unsigned int</type>
      <name>initializeDataOnPatchLevel</name>
      <anchorfile>class_i_b_t_k_1_1_l_init_strategy.html</anchorfile>
      <anchor>ab870b40b60125539b6477d6f00c0b9fb</anchor>
      <arglist>(int lag_node_index_idx, unsigned int global_index_offset, unsigned int local_index_offset, SAMRAI::tbox::Pointer&lt; LData &gt; X_data, SAMRAI::tbox::Pointer&lt; LData &gt; U_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, LDataManager *l_data_manager)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual unsigned int</type>
      <name>initializeMassDataOnPatchLevel</name>
      <anchorfile>class_i_b_t_k_1_1_l_init_strategy.html</anchorfile>
      <anchor>a64f8832d631c72151049651308bd3f98</anchor>
      <arglist>(unsigned int global_index_offset, unsigned int local_index_offset, SAMRAI::tbox::Pointer&lt; LData &gt; M_data, SAMRAI::tbox::Pointer&lt; LData &gt; K_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual unsigned int</type>
      <name>initializeDirectorDataOnPatchLevel</name>
      <anchorfile>class_i_b_t_k_1_1_l_init_strategy.html</anchorfile>
      <anchor>a3ce5fb1f9eeb14d6471f884a6e116104</anchor>
      <arglist>(unsigned int global_index_offset, unsigned int local_index_offset, SAMRAI::tbox::Pointer&lt; LData &gt; D_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, LDataManager *l_data_manager)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>tagCellsForInitialRefinement</name>
      <anchorfile>class_i_b_t_k_1_1_l_init_strategy.html</anchorfile>
      <anchor>a8a41a9e8ba9ba29131b5c239c65fc82d</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LMarker</name>
    <filename>class_i_b_t_k_1_1_l_marker.html</filename>
    <member kind="function">
      <type></type>
      <name>LMarker</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a735869610eb1bbd23c0ddae5a356620a</anchor>
      <arglist>(int idx=-1, const Point &amp;X=Point::Zero(), const Vector &amp;U=Vector::Zero(), const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_offset=SAMRAI::hier::IntVector&lt; NDIM &gt;(0))</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LMarker</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a54be15b897fca3169871fe9d23dd15e5</anchor>
      <arglist>(const LMarker &amp;from)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LMarker</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>ae44165ab7d8ec5695bdd3d474a00d2da</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LMarker</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a3d5e2ec77e26c3dae68096ea8c5b9a25</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>LMarker &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a65e3bc4f2bbcac9f3a88f601b7ef391b</anchor>
      <arglist>(const LMarker &amp;that)</arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>getIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a4dc54060ba6fd737e4117cba6e458473</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>getIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a2bfc8ec8d1769c53b6122afc3780bcc9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a78c6587e8a319e7ca2a4a54567eb9837</anchor>
      <arglist>(int idx)</arglist>
    </member>
    <member kind="function">
      <type>const Point &amp;</type>
      <name>getPosition</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>acdbc860b684ad37d0d93bf6d36eb0e3a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Point &amp;</type>
      <name>getPosition</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a08172f4f50c413efcb180c9f3233e782</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPosition</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>aa1b201819778ba5aaa48a97ae1832f0e</anchor>
      <arglist>(const Point &amp;X)</arglist>
    </member>
    <member kind="function">
      <type>const Vector &amp;</type>
      <name>getVelocity</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a30397970d17f70e96077e2022a335d4f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Vector &amp;</type>
      <name>getVelocity</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a01a6bad93fccc66a37fc69e580e88d0c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setVelocity</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a02b2cd3f087c2843ab4d4d59a6df46bb</anchor>
      <arglist>(const Vector &amp;U)</arglist>
    </member>
    <member kind="function">
      <type>const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getPeriodicOffset</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>ac43bf4b5df659e4e4126da06a371a853</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPeriodicOffset</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>ae60abf3ae724881131b816df59903aff</anchor>
      <arglist>(const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>copySourceItem</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a5d091c6c5d38584a2ade48961070bced</anchor>
      <arglist>(const SAMRAI::hier::Index&lt; NDIM &gt; &amp;src_index, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset, const LMarker &amp;src_item)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>afbc1718639897fa639463d862b250ddb</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>af5a8aa4ae7178652c71af8bb6f15aac3</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>unpackStream</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker.html</anchorfile>
      <anchor>a0162fdbfd897bd09650106eb04dc9a7f</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LMarkerCoarsen</name>
    <filename>class_i_b_t_k_1_1_l_marker_coarsen.html</filename>
    <base>CoarsenOperator&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>LMarkerCoarsen</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_coarsen.html</anchorfile>
      <anchor>ae549cf79965058e0559a7fed39b4a41e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LMarkerCoarsen</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_coarsen.html</anchorfile>
      <anchor>a559c7a31b28d30d187ac543c69fb5189</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>findCoarsenOperator</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_coarsen.html</anchorfile>
      <anchor>a49b05a0fa5fea079268dedf486c154c9</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;var, const std::string &amp;op_name) const </arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getOperatorName</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_coarsen.html</anchorfile>
      <anchor>ad5e042003aa635a46fdb0ec3a7a8781f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getOperatorPriority</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_coarsen.html</anchorfile>
      <anchor>a23176681a138a189534802d9d6d4aa08</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_coarsen.html</anchorfile>
      <anchor>a6b16a24a2c977df4be969da093934a29</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>coarsen</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_coarsen.html</anchorfile>
      <anchor>a7838f61ec7a034aa0f963095ea40a9b8</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, int dst_component, int src_component, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;coarse_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LMarkerRefine</name>
    <filename>class_i_b_t_k_1_1_l_marker_refine.html</filename>
    <base>RefineOperator&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>LMarkerRefine</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_refine.html</anchorfile>
      <anchor>ac193622a264a306552eb347234dba3cf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LMarkerRefine</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_refine.html</anchorfile>
      <anchor>a35f3d3032a5663eb26f5df802555bd04</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>findRefineOperator</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_refine.html</anchorfile>
      <anchor>ae42ceb31c70bbb4e3cd94a34a13d97a6</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;var, const std::string &amp;op_name) const </arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getOperatorName</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_refine.html</anchorfile>
      <anchor>af466d8eeb3e990701bd2c238c5a38e74</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getOperatorPriority</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_refine.html</anchorfile>
      <anchor>a33a6f27a2b5054f79fd35dc20bd6a1cf</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_refine.html</anchorfile>
      <anchor>adfe0620335b2720871d1780473f1d811</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>refine</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_refine.html</anchorfile>
      <anchor>aece6d96266a09a1a7535a269723ec7b0</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, int dst_component, int src_component, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LMarkerUtilities</name>
    <filename>class_i_b_t_k_1_1_l_marker_utilities.html</filename>
    <member kind="function" static="yes">
      <type>static unsigned int</type>
      <name>readMarkerPositions</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_utilities.html</anchorfile>
      <anchor>a4bf82c57b944895b67ad9f9c9fe21c48</anchor>
      <arglist>(std::vector&lt; Point &gt; &amp;mark_init_posns, const std::string &amp;mark_input_file_name, SAMRAI::tbox::Pointer&lt; SAMRAI::geom::CartesianGridGeometry&lt; NDIM &gt; &gt; grid_geom)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>eulerStep</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_utilities.html</anchorfile>
      <anchor>a6a3f1dade9fd6eaaa07b8790b6869e92</anchor>
      <arglist>(int mark_current_idx, int mark_new_idx, int u_current_idx, double dt, const std::string &amp;weighting_fcn, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>midpointStep</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_utilities.html</anchorfile>
      <anchor>a801fe887049217bbb8a269830ff97b2c</anchor>
      <arglist>(int mark_current_idx, int mark_new_idx, int u_half_idx, double dt, const std::string &amp;weighting_fcn, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>trapezoidalStep</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_utilities.html</anchorfile>
      <anchor>a7cd7d6724471d8d37f1c5b71d8bb1f1c</anchor>
      <arglist>(int mark_current_idx, int mark_new_idx, int u_new_idx, double dt, const std::string &amp;weighting_fcn, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>collectMarkersOnPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_utilities.html</anchorfile>
      <anchor>abd8f77fd440626bf51a9409cad2835e0</anchor>
      <arglist>(int mark_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>initializeMarkersOnLevel</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_utilities.html</anchorfile>
      <anchor>a67be6077b09db363aff6200abcf9bb0d</anchor>
      <arglist>(int mark_idx, const std::vector&lt; Point &gt; &amp;mark_init_posns, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>pruneInvalidMarkers</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_utilities.html</anchorfile>
      <anchor>a8e9f6c31f5548bad47ea8cc3a0899f02</anchor>
      <arglist>(int mark_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static unsigned int</type>
      <name>countMarkers</name>
      <anchorfile>class_i_b_t_k_1_1_l_marker_utilities.html</anchorfile>
      <anchor>a40d6c372e92580cf30bb33e269829168</anchor>
      <arglist>(int mark_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_ln=-1, int finest_ln=-1)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LMesh</name>
    <filename>class_i_b_t_k_1_1_l_mesh.html</filename>
    <member kind="function">
      <type></type>
      <name>LMesh</name>
      <anchorfile>class_i_b_t_k_1_1_l_mesh.html</anchorfile>
      <anchor>af7828fb655820bd5128f23700927bb33</anchor>
      <arglist>(const std::string &amp;object_name, const std::vector&lt; LNode * &gt; &amp;local_nodes, const std::vector&lt; LNode * &gt; &amp;ghost_nodes)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LMesh</name>
      <anchorfile>class_i_b_t_k_1_1_l_mesh.html</anchorfile>
      <anchor>abbd71895baa753b0a14522aa3507bb1c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; LNode * &gt; &amp;</type>
      <name>getLocalNodes</name>
      <anchorfile>class_i_b_t_k_1_1_l_mesh.html</anchorfile>
      <anchor>a0ecb9d65bf1ebbf63260dfd8954c0b09</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; LNode * &gt; &amp;</type>
      <name>getGhostNodes</name>
      <anchorfile>class_i_b_t_k_1_1_l_mesh.html</anchorfile>
      <anchor>a2a0b8b727e895ce451dff8d96abcfd3f</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LNode</name>
    <filename>class_i_b_t_k_1_1_l_node.html</filename>
    <base>IBTK::LNodeIndex</base>
    <member kind="function">
      <type></type>
      <name>LNode</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a67364b437e5aadaf6df2e1e0a3ae1e8d</anchor>
      <arglist>(int lagrangian_nidx=-1, int global_petsc_nidx=-1, int local_petsc_nidx=-1, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_offset=SAMRAI::hier::IntVector&lt; NDIM &gt;(0), const Vector &amp;periodic_displacement=Vector::Zero(), const std::vector&lt; SAMRAI::tbox::Pointer&lt; Streamable &gt; &gt; &amp;node_data=std::vector&lt; SAMRAI::tbox::Pointer&lt; Streamable &gt; &gt;())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LNode</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>aa5d5f886c30a266e0c69a07e54d8cac6</anchor>
      <arglist>(const LNode &amp;from)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LNode</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a816131e0bedc9b9dac7d9d03ecd74fb6</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LNode</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a28840b99b6c23150dad683fd8824ee38</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>LNode &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a9949b3243de30d8211d2f6de755a2551</anchor>
      <arglist>(const LNode &amp;that)</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; SAMRAI::tbox::Pointer&lt; Streamable &gt; &gt; &amp;</type>
      <name>getNodeData</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a06e3b7ed1be61f68ea30eb7a6e91cb78</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setNodeData</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a66b2d5805dcab2d3e98d5a1583589c76</anchor>
      <arglist>(const std::vector&lt; SAMRAI::tbox::Pointer&lt; Streamable &gt; &gt; &amp;node_data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>appendNodeDataItem</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>aab7bdbab4eda6c98fab21c8774fe37d1</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; Streamable &gt; &amp;node_data_item)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>removeNodeDataItem</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>aca29bf3d94f5828469504980ca463bc7</anchor>
      <arglist>(const SAMRAI::tbox::Pointer&lt; Streamable &gt; &amp;node_data_item)</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>getNodeDataItem</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a7ad31f89c6cb4914325436fabf355ecf</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; T * &gt;</type>
      <name>getNodeDataVector</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a823f68dd9b99896ca97cfc129b71c468</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerPeriodicShift</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a3941ca1cb32e91b2670cd63b19cbeec6</anchor>
      <arglist>(const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset, const Vector &amp;displacement)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>copySourceItem</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a663d78d9a33aeb205231d6a389fd1f68</anchor>
      <arglist>(const SAMRAI::hier::Index&lt; NDIM &gt; &amp;src_index, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset, const LNodeIndex &amp;src_item)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a9a71718b09f530b493a4a287a78243dd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>ad684817bf4792edc70c832d10767669a</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>unpackStream</name>
      <anchorfile>class_i_b_t_k_1_1_l_node.html</anchorfile>
      <anchor>a7f99573389a37303acbbab1618dd4e12</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LNodeIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>a43c046da44c17343c2209d81883ea1c1</anchor>
      <arglist>(int lagrangian_nidx=-1, int global_petsc_nidx=-1, int local_petsc_nidx=-1, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic_offset=SAMRAI::hier::IntVector&lt; NDIM &gt;(0), const Vector &amp;periodic_displacement=Vector::Zero())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LNodeIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>ada0369c92681faf38c2ed4b70d078ee4</anchor>
      <arglist>(const LNodeIndex &amp;from)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LNodeIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>abf0dda362ea223c7dd4542cb194b76b5</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LNodeIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>aab629f95dc4a4e93b488c4af6157e3bb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>LNodeIndex &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>a2d5bb428a36965451e092b9c86f066fc</anchor>
      <arglist>(const LNodeIndex &amp;that)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getLagrangianIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>af32ece0cfabd083e59865a7814a61e6b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setLagrangianIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>ae15f0387a85dd2214ab854f6e1ae8e91</anchor>
      <arglist>(int lagrangian_nidx)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getGlobalPETScIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>aaa5bf6c4e1824eb42e65e8ff92b0ef1e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setGlobalPETScIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>a7aa0995551e73282d9cfa8d623a6506f</anchor>
      <arglist>(int global_petsc_nidx)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getLocalPETScIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>abbe4b4b40434d137379330463f42609f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setLocalPETScIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>aec9cce8981d322033597386d99a34060</anchor>
      <arglist>(int local_petsc_nidx)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getPeriodicOffset</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>ac9dede00d2a403446f8b7c9682966187</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const Vector &amp;</type>
      <name>getPeriodicDisplacement</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>a90d7a52671133b1bf7f95e113c54ba9a</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LNodeIndex</name>
    <filename>class_i_b_t_k_1_1_l_node_index.html</filename>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerPeriodicShift</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>abaedc2c9c88fb55f341d62303f1f1b30</anchor>
      <arglist>(const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset, const Vector &amp;displacement)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>copySourceItem</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>a55b89fcef616cce2e04c9c4be2b4acb7</anchor>
      <arglist>(const SAMRAI::hier::Index&lt; NDIM &gt; &amp;src_index, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset, const LNodeIndex &amp;src_item)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>ab36b1ff8fab319f19395883c7ac3f8f8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>a1274b7cd3f4ae677a3be4aa3b9c72b02</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>unpackStream</name>
      <anchorfile>class_i_b_t_k_1_1_l_node_index.html</anchorfile>
      <anchor>a2bbe9a9421566a4b7460e699be32d3f2</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>IBTK::LNodeIndexGlobalPETScIndexComp</name>
    <filename>struct_i_b_t_k_1_1_l_node_index_global_p_e_t_sc_index_comp.html</filename>
    <base>binary_function&lt; const LNodeIndex &amp;, const LNodeIndex &amp;, bool &gt;</base>
    <base>binary_function&lt; const LNodeIndex *, const LNodeIndex *, bool &gt;</base>
  </compound>
  <compound kind="struct">
    <name>IBTK::LNodeIndexGlobalPETScIndexEqual</name>
    <filename>struct_i_b_t_k_1_1_l_node_index_global_p_e_t_sc_index_equal.html</filename>
    <base>binary_function&lt; const LNodeIndex &amp;, const LNodeIndex &amp;, bool &gt;</base>
    <base>binary_function&lt; const LNodeIndex *, const LNodeIndex *, bool &gt;</base>
  </compound>
  <compound kind="struct">
    <name>IBTK::LNodeIndexLagrangianIndexComp</name>
    <filename>struct_i_b_t_k_1_1_l_node_index_lagrangian_index_comp.html</filename>
    <base>binary_function&lt; const LNodeIndex &amp;, const LNodeIndex &amp;, bool &gt;</base>
    <base>binary_function&lt; const LNodeIndex *, const LNodeIndex *, bool &gt;</base>
  </compound>
  <compound kind="struct">
    <name>IBTK::LNodeIndexLagrangianIndexEqual</name>
    <filename>struct_i_b_t_k_1_1_l_node_index_lagrangian_index_equal.html</filename>
    <base>binary_function&lt; const LNodeIndex &amp;, const LNodeIndex &amp;, bool &gt;</base>
    <base>binary_function&lt; const LNodeIndex *, const LNodeIndex *, bool &gt;</base>
  </compound>
  <compound kind="struct">
    <name>IBTK::LNodeIndexLocalPETScIndexComp</name>
    <filename>struct_i_b_t_k_1_1_l_node_index_local_p_e_t_sc_index_comp.html</filename>
    <base>binary_function&lt; const LNodeIndex &amp;, const LNodeIndex &amp;, bool &gt;</base>
    <base>binary_function&lt; const LNodeIndex *, const LNodeIndex *, bool &gt;</base>
  </compound>
  <compound kind="struct">
    <name>IBTK::LNodeIndexLocalPETScIndexEqual</name>
    <filename>struct_i_b_t_k_1_1_l_node_index_local_p_e_t_sc_index_equal.html</filename>
    <base>binary_function&lt; const LNodeIndex &amp;, const LNodeIndex &amp;, bool &gt;</base>
    <base>binary_function&lt; const LNodeIndex *, const LNodeIndex *, bool &gt;</base>
  </compound>
  <compound kind="class">
    <name>IBTK::LNodeIndexPosnComp</name>
    <filename>class_i_b_t_k_1_1_l_node_index_posn_comp.html</filename>
    <base protection="private">binary_function&lt; const LNodeIndex &amp;, const LNodeIndex &amp;, bool &gt;</base>
    <base protection="private">binary_function&lt; const LNodeIndex *, const LNodeIndex *, bool &gt;</base>
  </compound>
  <compound kind="class">
    <name>IBTK::LNodeIndexPosnEqual</name>
    <filename>class_i_b_t_k_1_1_l_node_index_posn_equal.html</filename>
    <base protection="private">binary_function&lt; const LNodeIndex &amp;, const LNodeIndex &amp;, bool &gt;</base>
    <base protection="private">binary_function&lt; const LNodeIndex *, const LNodeIndex *, bool &gt;</base>
  </compound>
  <compound kind="class">
    <name>IBTK::LSet</name>
    <filename>class_i_b_t_k_1_1_l_set.html</filename>
    <templarg></templarg>
    <member kind="typedef">
      <type>std::vector&lt; SAMRAI::tbox::Pointer&lt; T &gt; &gt;</type>
      <name>DataSet</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a3802a4a6d519581c84ddc3225a2d244f</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DataSet::value_type</type>
      <name>value_type</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a1ba3ed5718695de898194c0d493d0425</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DataSet::pointer</type>
      <name>pointer</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a07142efc6444df7055a65737cabc4d59</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DataSet::reference</type>
      <name>reference</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a9db66579f8036cc855ca8da769a01fed</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DataSet::const_reference</type>
      <name>const_reference</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a5d4cfab5f706fec7f3cd0684796ae084</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DataSet::size_type</type>
      <name>size_type</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a4eb6a3fc4982b3cd1a505d974585f089</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DataSet::difference_type</type>
      <name>difference_type</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>ad885856d5847e40d127a065f09e68b81</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DataSet::iterator</type>
      <name>iterator</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a0265d09ba455081a0c6b55001fdc207e</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>DataSet::const_iterator</type>
      <name>const_iterator</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a9ba3db78eecd78379a3ab8982b86e65d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LSet</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a4146f5be948079877b1366eb04449095</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LSet</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a415a27fa2203b964b047073f7ab3e98c</anchor>
      <arglist>(const LSet &amp;from)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LSet</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a9a68ded380148d10ed758deb79dae1f5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>LSet &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a6083c9c6c3a88660574c0b3f024e2588</anchor>
      <arglist>(const LSet &amp;that)</arglist>
    </member>
    <member kind="function">
      <type>reference</type>
      <name>operator[]</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>af8ae24a4fe17d5c52600ad586867ca2c</anchor>
      <arglist>(size_type n)</arglist>
    </member>
    <member kind="function">
      <type>const_reference</type>
      <name>operator[]</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a8a36bda1063e41b72d64022f06c6df0c</anchor>
      <arglist>(size_type n) const </arglist>
    </member>
    <member kind="function">
      <type>const_iterator</type>
      <name>begin</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>aea9013b0f68b25d65a636e79f70855f5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>iterator</type>
      <name>begin</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>ae5b04cafe29181596d60735db7f27304</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const_iterator</type>
      <name>end</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a187b1de94f7f3988f20122b605f396ea</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>iterator</type>
      <name>end</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a6912a27ad26bdccd0989686d3304ee12</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>size_type</type>
      <name>size</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a6830b0df6ea081165fd1969e61c47ee8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>empty</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a3aef386ed2c8c0eb881677b58133855b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>push_back</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a21e787af111b9fdacfa1068d6bb5a8c2</anchor>
      <arglist>(const value_type &amp;value)</arglist>
    </member>
    <member kind="function">
      <type>iterator</type>
      <name>insert</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a9c4c8c8ef8548459b40d1ddc73b6e503</anchor>
      <arglist>(iterator pos, const typename LSet&lt; T &gt;::value_type &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>insert</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>acbb21370602bed1115753413fa48aa2a</anchor>
      <arglist>(iterator pos, InputIterator first, InputIterator last)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>insert</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a479655f0de2c9a06058a1b702fdea006</anchor>
      <arglist>(iterator pos, size_type n, const typename LSet&lt; T &gt;::value_type &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>const DataSet &amp;</type>
      <name>getDataSet</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>aa5404638389f06a3eeadc9e58dfe570c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DataSet &amp;</type>
      <name>getDataSet</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a5c1d4e185eb7113609545c2d92d4b9a7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataSet</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a477b80d3f75a2ea243c3837d36de6bba</anchor>
      <arglist>(const DataSet &amp;set)</arglist>
    </member>
    <member kind="function">
      <type>const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getPeriodicOffset</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>afdd77bdff249173dae98fd2bd4f3d1b5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPeriodicOffset</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a9428c005c09c5676343e1dbb8468085d</anchor>
      <arglist>(const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>copySourceItem</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a105e9f4d3b8f0184d16a998254c1318d</anchor>
      <arglist>(const SAMRAI::hier::Index&lt; NDIM &gt; &amp;src_index, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset, const LSet &amp;src_item)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a6d1a2c22206f3d02b460eafc923be93a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a701fd94e0bde013548b03026c62af6fc</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unpackStream</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>ad2cb4c07dc0e694813885d1f29047942</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>abfeb053e80e8626d1ab0fbde8ab3620a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; database)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>getFromDatabase</name>
      <anchorfile>class_i_b_t_k_1_1_l_set.html</anchorfile>
      <anchor>a296661c2a469fef19d671c99181ca765</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; database)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LSetData</name>
    <filename>class_i_b_t_k_1_1_l_set_data.html</filename>
    <templarg></templarg>
    <base>IndexData&lt; NDIM, LSet&lt; T &gt;, SAMRAI::pdat::CellGeometry&lt; NDIM &gt; &gt;</base>
  </compound>
  <compound kind="class">
    <name>IBTK::LSetDataFactory</name>
    <filename>class_i_b_t_k_1_1_l_set_data_factory.html</filename>
    <templarg></templarg>
    <base>IndexDataFactory&lt; NDIM, LSet&lt; T &gt;, SAMRAI::pdat::CellGeometry&lt; NDIM &gt; &gt;</base>
  </compound>
  <compound kind="class">
    <name>IBTK::LSetDataIterator</name>
    <filename>class_i_b_t_k_1_1_l_set_data_iterator.html</filename>
    <templarg></templarg>
    <member kind="function">
      <type></type>
      <name>LSetDataIterator</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_iterator.html</anchorfile>
      <anchor>a3952f8c02b2bcac24eda1d4759967eee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LSetDataIterator</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_iterator.html</anchorfile>
      <anchor>a903f8815fb726741493eb82adf98bdcb</anchor>
      <arglist>(const LSetDataIterator &amp;that)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LSetDataIterator</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_iterator.html</anchorfile>
      <anchor>a9fa06cf18f03ac95e6b632f133715a4d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>LSetDataIterator&lt; T &gt; &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_iterator.html</anchorfile>
      <anchor>af6d38cb35db53346dd422d01ba15c2e5</anchor>
      <arglist>(const LSetDataIterator&lt; T &gt; &amp;that)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_iterator.html</anchorfile>
      <anchor>ad5d9a3d021ca293f9aa948b8e5c72753</anchor>
      <arglist>(const LSetDataIterator&lt; T &gt; &amp;that)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator!=</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_iterator.html</anchorfile>
      <anchor>a8a9b8abda24af3d03409a41de3a78dfc</anchor>
      <arglist>(const LSetDataIterator&lt; T &gt; &amp;that)</arglist>
    </member>
    <member kind="function">
      <type>LSetDataIterator&lt; T &gt; &amp;</type>
      <name>operator++</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_iterator.html</anchorfile>
      <anchor>a6de5755c2f7fdf40037f9bf93bb71adb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>LSetDataIterator&lt; T &gt;</type>
      <name>operator++</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_iterator.html</anchorfile>
      <anchor>abe8c9c6570509aef0bb9866ca67f134f</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>LSet&lt; T &gt;::value_type &amp;</type>
      <name>operator*</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_iterator.html</anchorfile>
      <anchor>aa1ef3536f0946934d82da3ffa37f8a5a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>LSet&lt; T &gt;::value_type &amp;</type>
      <name>getDataItem</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_iterator.html</anchorfile>
      <anchor>ace71f2481e6805a994e5546c061a8d6e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const SAMRAI::hier::Index&lt; NDIM &gt; &amp;</type>
      <name>getCellIndex</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_data_iterator.html</anchorfile>
      <anchor>a65e9ab27910fdd187e096882a502a281</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LSetVariable</name>
    <filename>class_i_b_t_k_1_1_l_set_variable.html</filename>
    <templarg></templarg>
    <base>Variable&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>LSetVariable</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_variable.html</anchorfile>
      <anchor>a13b2e4c3bbafa97fff0d67a6304913e3</anchor>
      <arglist>(const std::string &amp;name)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LSetVariable</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_variable.html</anchorfile>
      <anchor>adb302bbde3c70fd4ec0fc25e26a2406d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dataLivesOnPatchBorder</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_variable.html</anchorfile>
      <anchor>a7501180ba3a38db4253dcab7d416f282</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>fineBoundaryRepresentsVariable</name>
      <anchorfile>class_i_b_t_k_1_1_l_set_variable.html</anchorfile>
      <anchor>af57f4ee9c748d46853ec3d32e99f4ee3</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LSiloDataWriter</name>
    <filename>class_i_b_t_k_1_1_l_silo_data_writer.html</filename>
    <base>SAMRAI::tbox::Serializable</base>
    <member kind="function">
      <type></type>
      <name>LSiloDataWriter</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>a8c5fea0642d57620ae18fcbb558ba7d7</anchor>
      <arglist>(const std::string &amp;object_name, const std::string &amp;dump_directory_name, bool register_for_restart=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LSiloDataWriter</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>aa1e86ff2d4ab8c83f8bc5a405a86d787</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerMarkerCloud</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>a4e67a27464397ee80b97a59f6e6a7a5f</anchor>
      <arglist>(const std::string &amp;name, int nmarks, int first_lag_idx, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLogicallyCartesianBlock</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>ae1a633968cc34df44449f852909c70e7</anchor>
      <arglist>(const std::string &amp;name, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;nelem, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;periodic, int first_lag_idx, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLogicallyCartesianMultiblock</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>afe63dda494dd820afc7db23bd87234e1</anchor>
      <arglist>(const std::string &amp;name, const std::vector&lt; SAMRAI::hier::IntVector&lt; NDIM &gt; &gt; &amp;nelem, const std::vector&lt; SAMRAI::hier::IntVector&lt; NDIM &gt; &gt; &amp;periodic, const std::vector&lt; int &gt; &amp;first_lag_idx, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerUnstructuredMesh</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>ab7f2705bcc83fe95a1c76b40b9d22b59</anchor>
      <arglist>(const std::string &amp;name, const std::multimap&lt; int, std::pair&lt; int, int &gt; &gt; &amp;edge_map, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerCoordsData</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>a50ba8067a423c4268a59a040f63b2c60</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LData &gt; coords_data, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerVariableData</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>a1e4c72b046b963a7ae9c02e0a0e22f73</anchor>
      <arglist>(const std::string &amp;var_name, SAMRAI::tbox::Pointer&lt; LData &gt; var_data, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerVariableData</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>adc0456f7f20ba766211ead4de6ae69ca</anchor>
      <arglist>(const std::string &amp;var_name, SAMRAI::tbox::Pointer&lt; LData &gt; var_data, int start_depth, int var_depth, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLagrangianAO</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>a98947c98e7838127a6a75c5b168f50b7</anchor>
      <arglist>(AO &amp;ao, int level_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerLagrangianAO</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>a2dd683a7d04e34fc85dd683b511cb1ce</anchor>
      <arglist>(std::vector&lt; AO &gt; &amp;ao, int coarsest_ln, int finest_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>writePlotData</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>a5d4aab895aec4c22c12017638e3ed26b</anchor>
      <arglist>(int time_step_number, double simulation_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putToDatabase</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>a6903a7bae3ac65037d77cd9a93bdd4ff</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; db)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPatchHierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>ac1bee41fe4c3ab592dd88db2758256af</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetLevels</name>
      <anchorfile>class_i_b_t_k_1_1_l_silo_data_writer.html</anchorfile>
      <anchor>a7c618d3bed71024e561344e1e1347da3</anchor>
      <arglist>(int coarsest_ln, int finest_ln)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LTransaction</name>
    <filename>class_i_b_t_k_1_1_l_transaction.html</filename>
    <templarg></templarg>
    <base>SAMRAI::tbox::Transaction</base>
    <class kind="class">IBTK::LTransaction::LTransactionComponent</class>
    <member kind="function">
      <type></type>
      <name>LTransaction</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>a906ce39af7098cead918baf9a1ac388c</anchor>
      <arglist>(int src_proc, int dst_proc)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LTransaction</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>af2d4d0b638e4a5f666e14fe10409fe15</anchor>
      <arglist>(int src_proc, int dst_proc, const std::vector&lt; LTransactionComponent &gt; &amp;src_item_set)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LTransaction</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>ab163e82a4fbd41cf2efa683d7d31a2a3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; LTransactionComponent &gt; &amp;</type>
      <name>getSourceData</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>afeabd8f74de498e6b336c7f0408c2c49</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; LTransactionComponent &gt; &amp;</type>
      <name>getDestinationData</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>a8644886c4da91812d2096fee7e674041</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>canEstimateIncomingMessageSize</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>a02cbf7b1ba19513eb13cc20e61b84791</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>computeIncomingMessageSize</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>a773e4bc9845f340f396dc3d29b1ae5b8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>computeOutgoingMessageSize</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>a4acc4d8f65a40b0f7c0f7cb8da155553</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getSourceProcessor</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>ae487a226487800ec2d98ce30b10923ae</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getDestinationProcessor</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>a688a4f8348e1879242eaf71889f999e6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>a6375e94159179bef51fc9375617ce409</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>unpackStream</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>abc6bbf8267573d4c7bfcf9c69a5e6fb6</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>copyLocalData</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>a49b41e9c3dd9268d94857a763927fd7b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>printClassData</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction.html</anchorfile>
      <anchor>a89e0254d0c4c7ae5d0920db4a32a79e0</anchor>
      <arglist>(std::ostream &amp;stream) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::LTransaction::LTransactionComponent</name>
    <filename>class_i_b_t_k_1_1_l_transaction_1_1_l_transaction_component.html</filename>
    <member kind="function">
      <type></type>
      <name>LTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction_1_1_l_transaction_component.html</anchorfile>
      <anchor>aa03dc363a356ffac22e6a693fd896d03</anchor>
      <arglist>(const typename LSet&lt; T &gt;::value_type &amp;item=NULL, const Point &amp;posn=Point::Zero())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction_1_1_l_transaction_component.html</anchorfile>
      <anchor>a5e5583df1738685435aa47e7505fb4da</anchor>
      <arglist>(const LTransactionComponent &amp;from)</arglist>
    </member>
    <member kind="function">
      <type>LTransactionComponent &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction_1_1_l_transaction_component.html</anchorfile>
      <anchor>a052f1f81d97635c12479c0d82ee812a7</anchor>
      <arglist>(const LTransactionComponent &amp;that)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_l_transaction_1_1_l_transaction_component.html</anchorfile>
      <anchor>ab90613bd124e290465c9cc5fb74907d1</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::muParserCartGridFunction</name>
    <filename>class_i_b_t_k_1_1mu_parser_cart_grid_function.html</filename>
    <base>IBTK::CartGridFunction</base>
    <member kind="function">
      <type></type>
      <name>muParserCartGridFunction</name>
      <anchorfile>class_i_b_t_k_1_1mu_parser_cart_grid_function.html</anchorfile>
      <anchor>a2eaafa0de9f6f9d8a6155efe0370a503</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, SAMRAI::tbox::Pointer&lt; SAMRAI::geom::CartesianGridGeometry&lt; NDIM &gt; &gt; grid_geom)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~muParserCartGridFunction</name>
      <anchorfile>class_i_b_t_k_1_1mu_parser_cart_grid_function.html</anchorfile>
      <anchor>aa362b58dbf2df6a413e264d2789316a1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isTimeDependent</name>
      <anchorfile>class_i_b_t_k_1_1mu_parser_cart_grid_function.html</anchorfile>
      <anchor>adccb1f27928888b24dcf7ec3a4b02c17</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setDataOnPatch</name>
      <anchorfile>class_i_b_t_k_1_1mu_parser_cart_grid_function.html</anchorfile>
      <anchor>a26b2211d0023632ba066a947b097dc10</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; var, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, double data_time, bool initial_time=false, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt;(NULL))</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::muParserRobinBcCoefs</name>
    <filename>class_i_b_t_k_1_1mu_parser_robin_bc_coefs.html</filename>
    <base>RobinBcCoefStrategy&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>muParserRobinBcCoefs</name>
      <anchorfile>class_i_b_t_k_1_1mu_parser_robin_bc_coefs.html</anchorfile>
      <anchor>abc499e2f29b42bfc85e0cb6c17d466e2</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, SAMRAI::tbox::Pointer&lt; SAMRAI::geom::CartesianGridGeometry&lt; NDIM &gt; &gt; grid_geom)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~muParserRobinBcCoefs</name>
      <anchorfile>class_i_b_t_k_1_1mu_parser_robin_bc_coefs.html</anchorfile>
      <anchor>ae0b4d0f263bb87898a18c4ac3be5766b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setBcCoefs</name>
      <anchorfile>class_i_b_t_k_1_1mu_parser_robin_bc_coefs.html</anchorfile>
      <anchor>a765871401da1ecfa1abfacf1d2da7c57</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;acoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;bcoef_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::ArrayData&lt; NDIM, double &gt; &gt; &amp;gcoef_data, const SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Variable&lt; NDIM &gt; &gt; &amp;variable, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, const SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &amp;bdry_box, double fill_time=0.0) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::NewtonKrylovSolver</name>
    <filename>class_i_b_t_k_1_1_newton_krylov_solver.html</filename>
    <base>IBTK::GeneralSolver</base>
    <member kind="function">
      <type></type>
      <name>NewtonKrylovSolver</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>a43c8d1903b5807604f23c2b475f09dd7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~NewtonKrylovSolver</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>a099eb2c89a8681ae93e4144f07729254</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHierarchyMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>ad69c8d0a8e32143284273bc44fe4f870</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; HierarchyMathOps &gt; hier_math_ops)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setHomogeneousBc</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>ae1cf5d41fd9b64ee1037798c6242e807</anchor>
      <arglist>(bool homogeneous_bc)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSolutionTime</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>a49b3f9155f320d847cd547a287c6c0fe</anchor>
      <arglist>(double solution_time)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTimeInterval</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>ae87e5de6308f31758f2fdbfa002525d6</anchor>
      <arglist>(double current_time, double new_time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setOperator</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>af0821917f74cb297b1169df584a28432</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; GeneralOperator &gt; op)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual SAMRAI::tbox::Pointer&lt; GeneralOperator &gt;</type>
      <name>getOperator</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>a33a9320ae47929037cb47ecc89b3e0f7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt;</type>
      <name>getSolutionVector</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>af87e3c77baa4278c0fd908d8ed6b31d1</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt;</type>
      <name>getFunctionVector</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>a05fef714fa15c384c3e5c8f17e8c3e76</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setJacobian</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>a6cc3d4d5fa85b48372496c08af9946f9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; JacobianOperator &gt; J)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual SAMRAI::tbox::Pointer&lt; JacobianOperator &gt;</type>
      <name>getJacobian</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>a47b8ef65b040de5990cbce279c1b917b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual SAMRAI::tbox::Pointer&lt; KrylovLinearSolver &gt;</type>
      <name>getLinearSolver</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>a49cec8caded7ce116ad151d3f12e7bc2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setMaxEvaluations</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>ab39f3ce8007de4f34f17d7e9d980d1f1</anchor>
      <arglist>(int max_evaluations)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getMaxEvaluations</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>a5e29cf4ee819b64b41cf65d98d47dee5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setSolutionTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>aff08c3fb7991cae46ef29b8507d652a6</anchor>
      <arglist>(double solution_tol)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>getSolutionTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>a3c97a33743518f2d59ce51be86448eda</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getNumLinearIterations</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver.html</anchorfile>
      <anchor>a20fc1c95eea851a6df7403548dd8fcf1</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::NewtonKrylovSolverManager</name>
    <filename>class_i_b_t_k_1_1_newton_krylov_solver_manager.html</filename>
    <member kind="typedef">
      <type>SAMRAI::tbox::Pointer&lt; NewtonKrylovSolver &gt;(*</type>
      <name>SolverMaker</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver_manager.html</anchorfile>
      <anchor>a47111264ddf3f990ec869f6b2853ea66</anchor>
      <arglist>)(const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; NewtonKrylovSolver &gt;</type>
      <name>allocateSolver</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver_manager.html</anchorfile>
      <anchor>a2ca6b14fb8fb7117e5a4d55019fcba06</anchor>
      <arglist>(const std::string &amp;solver_type, const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerSolverFactoryFunction</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver_manager.html</anchorfile>
      <anchor>a1b230157273c9b82759ac6d2a33f265c</anchor>
      <arglist>(const std::string &amp;solver_type, SolverMaker solver_maker)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static NewtonKrylovSolverManager *</type>
      <name>getManager</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver_manager.html</anchorfile>
      <anchor>af12259abf6852de7aeba8b6910a134a5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>freeManager</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver_manager.html</anchorfile>
      <anchor>a710ca9df4737dbe9b3f835908a143ade</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>UNDEFINED</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver_manager.html</anchorfile>
      <anchor>aca1e7e88e818946a6f086eab97bc8f71</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver_manager.html</anchorfile>
      <anchor>addf8bbce547cf87789a5f4a488f23a65</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>NewtonKrylovSolverManager</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver_manager.html</anchorfile>
      <anchor>ad9dd2c701125d20d3575f0a84d03cbac</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>~NewtonKrylovSolverManager</name>
      <anchorfile>class_i_b_t_k_1_1_newton_krylov_solver_manager.html</anchorfile>
      <anchor>a6cf35bb1f8d657591815c7cdc13ab252</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::NodeDataSynchronization</name>
    <filename>class_i_b_t_k_1_1_node_data_synchronization.html</filename>
    <class kind="class">IBTK::NodeDataSynchronization::SynchronizationTransactionComponent</class>
    <member kind="function">
      <type></type>
      <name>NodeDataSynchronization</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization.html</anchorfile>
      <anchor>a1262097916ea4d2382d94d3ceca0d8e6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~NodeDataSynchronization</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization.html</anchorfile>
      <anchor>ae59a75e509fe83e6f2a339e6f8c2041f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization.html</anchorfile>
      <anchor>a829f88bdf1fa9647f2f7bc416732b7c9</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;transaction_comp, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization.html</anchorfile>
      <anchor>a0a96da5fa5343a80998b94dd077b21bb</anchor>
      <arglist>(const std::vector&lt; SynchronizationTransactionComponent &gt; &amp;transaction_comps, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization.html</anchorfile>
      <anchor>ae401ac915ffe4a4b17d775ff2741aaa1</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;transaction_comps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTransactionComponents</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization.html</anchorfile>
      <anchor>ac01c01c1ae2a6f2a77e93c46343dbc3b</anchor>
      <arglist>(const std::vector&lt; SynchronizationTransactionComponent &gt; &amp;transaction_comps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization.html</anchorfile>
      <anchor>a94a19d387fa3f73e531535ca351988f4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>synchronizeData</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization.html</anchorfile>
      <anchor>a8b5c13f8dcae17b1a027ef0aaf64e251</anchor>
      <arglist>(double fill_time)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::NodeDataSynchronization::SynchronizationTransactionComponent</name>
    <filename>class_i_b_t_k_1_1_node_data_synchronization_1_1_synchronization_transaction_component.html</filename>
    <member kind="function">
      <type></type>
      <name>SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>a8105ccc8f3f79b00b2d731fcb1eaadb9</anchor>
      <arglist>(int data_idx=-1, const std::string &amp;coarsen_op_name=&quot;NONE&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>ad3fc9deb764e35bf50f497378bb91b18</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;from)</arglist>
    </member>
    <member kind="function">
      <type>SynchronizationTransactionComponent &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>a8eab296c59d47553bc4d8bb9f76c1f45</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;that)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_node_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>a2bff69ff47b02ab4f4b3f10c8febfde9</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::NodeSynchCopyFillPattern</name>
    <filename>class_i_b_t_k_1_1_node_synch_copy_fill_pattern.html</filename>
    <base>VariableFillPattern&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>NodeSynchCopyFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_node_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a6203ab066e2d31c4bdab962aa4bc33b4</anchor>
      <arglist>(unsigned int axis)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~NodeSynchCopyFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_node_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a0cf04ae655a27d95835c918f8a9e7f73</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BoxOverlap&lt; NDIM &gt; &gt;</type>
      <name>calculateOverlap</name>
      <anchorfile>class_i_b_t_k_1_1_node_synch_copy_fill_pattern.html</anchorfile>
      <anchor>ae90b3670f6b44ebeae4f5e68b2ed5546</anchor>
      <arglist>(const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;dst_geometry, const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;src_geometry, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;dst_patch_box, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;src_mask, bool overwrite_interior, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_node_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a0db3756492dd275d3116cfe29f0fdb72</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getPatternName</name>
      <anchorfile>class_i_b_t_k_1_1_node_synch_copy_fill_pattern.html</anchorfile>
      <anchor>af2e4e841f62698a9e424f6e3abd4be95</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::NormOps</name>
    <filename>class_i_b_t_k_1_1_norm_ops.html</filename>
    <member kind="function" static="yes">
      <type>static double</type>
      <name>L1Norm</name>
      <anchorfile>class_i_b_t_k_1_1_norm_ops.html</anchorfile>
      <anchor>a0af692737e1620bf0091dc4c45e8f46b</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; *samrai_vector, bool local_only=false)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static double</type>
      <name>L2Norm</name>
      <anchorfile>class_i_b_t_k_1_1_norm_ops.html</anchorfile>
      <anchor>a600796e7671d6809c9e1dacfce9bb0dc</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; *samrai_vector, bool local_only=false)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static double</type>
      <name>maxNorm</name>
      <anchorfile>class_i_b_t_k_1_1_norm_ops.html</anchorfile>
      <anchor>ada4b025a9bb501237458b739d7a59c16</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; *samrai_vector, bool local_only=false)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::ParallelEdgeMap</name>
    <filename>class_i_b_t_k_1_1_parallel_edge_map.html</filename>
    <member kind="function">
      <type></type>
      <name>ParallelEdgeMap</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_edge_map.html</anchorfile>
      <anchor>a44f9f2fbc83e4774ae4617067743d801</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ParallelEdgeMap</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_edge_map.html</anchorfile>
      <anchor>a1d8c3ef8eaf0d572fcd1d33884bd528f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>addEdge</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_edge_map.html</anchorfile>
      <anchor>acfab14e9e5fccc8745338c8dc3e2e885</anchor>
      <arglist>(const std::pair&lt; int, int &gt; &amp;link, int mastr_idx=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>removeEdge</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_edge_map.html</anchorfile>
      <anchor>a7801fe8f301ede8c258047995a34dd09</anchor>
      <arglist>(const std::pair&lt; int, int &gt; &amp;link, int mastr_idx=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>communicateData</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_edge_map.html</anchorfile>
      <anchor>aae4b3720740816d19dac15bd50bf1053</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::multimap&lt; int, std::pair&lt; int, int &gt; &gt; &amp;</type>
      <name>getEdgeMap</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_edge_map.html</anchorfile>
      <anchor>a8771fa89900f03fe8e01273c6a513338</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::ParallelMap</name>
    <filename>class_i_b_t_k_1_1_parallel_map.html</filename>
    <member kind="function">
      <type></type>
      <name>ParallelMap</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_map.html</anchorfile>
      <anchor>a13061aa96ed1b04ed51d39b81217f717</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ParallelMap</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_map.html</anchorfile>
      <anchor>a666956f8c51df91db5d0a7d2887ae59b</anchor>
      <arglist>(const ParallelMap &amp;from)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ParallelMap</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_map.html</anchorfile>
      <anchor>a889cfa924846eeab8b409ab401132ec1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>ParallelMap &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_map.html</anchorfile>
      <anchor>a648cb58357c2afe365070e6aee209d2a</anchor>
      <arglist>(const ParallelMap &amp;that)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>addItem</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_map.html</anchorfile>
      <anchor>af35855c13bdd9725d3cf14f0781b75f1</anchor>
      <arglist>(int key, SAMRAI::tbox::Pointer&lt; Streamable &gt; item)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>removeItem</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_map.html</anchorfile>
      <anchor>aa7cedad0d97be5991c37f42b4086aa45</anchor>
      <arglist>(int key)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>communicateData</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_map.html</anchorfile>
      <anchor>a719bf4a303a2e234826bf82428ae89c2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::map&lt; int, SAMRAI::tbox::Pointer&lt; Streamable &gt; &gt; &amp;</type>
      <name>getMap</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_map.html</anchorfile>
      <anchor>a6edf369ec3d56933376d82acf37204b7</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::ParallelSet</name>
    <filename>class_i_b_t_k_1_1_parallel_set.html</filename>
    <member kind="function">
      <type></type>
      <name>ParallelSet</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_set.html</anchorfile>
      <anchor>a88bb075b9a28b4cdc72f32a97061a972</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ParallelSet</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_set.html</anchorfile>
      <anchor>a5930f6d9139be8815ab788b7fd6799f7</anchor>
      <arglist>(const ParallelSet &amp;from)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ParallelSet</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_set.html</anchorfile>
      <anchor>a830748674961054a5a4510facb0540e9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>ParallelSet &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_set.html</anchorfile>
      <anchor>a7745c3fe93e624e961c337ab2e73aee9</anchor>
      <arglist>(const ParallelSet &amp;that)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>addItem</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_set.html</anchorfile>
      <anchor>a8b93801fc51d92c4274fb32aaa39e608</anchor>
      <arglist>(int key)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>removeItem</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_set.html</anchorfile>
      <anchor>a28cd4912f381f36349f1b369b7a6320a</anchor>
      <arglist>(int key)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>communicateData</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_set.html</anchorfile>
      <anchor>a13dbbd477aabd8c803589905b2eb01e4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::set&lt; int &gt; &amp;</type>
      <name>getSet</name>
      <anchorfile>class_i_b_t_k_1_1_parallel_set.html</anchorfile>
      <anchor>a90b0d096b10b7999757a62bd7605fe3c</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PatchMathOps</name>
    <filename>class_i_b_t_k_1_1_patch_math_ops.html</filename>
    <member kind="function">
      <type></type>
      <name>PatchMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>af18ed8639c1b2ea89f868f7e139226a1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PatchMathOps</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a1dc32cbf5d1301f8f14ad4efb619b79f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>aff8658ef85a6d6cc4ebc2e00535c0306</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ae542be793f11c94b37bbbba0c7f9db25</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a54690aa5faf4daa7cc23644d2070cd06</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ac36a3328f9dd1da26fd1871edf7ac251</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a4e24fcd7c3ae0f6c62124afdcaa5c606</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a35bcc1e885fcf437fbd3867516c0ace3</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>curl</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>abf6fc5d6f0e7f26cf1c8ef14ba6c7da9</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rot</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ac6d2eb74a341ccac4687a0faaac92b1e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rot</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ab8c8f14f4a2e9c3c49dafc58bdc5ce7a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rot</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ad8f1586f098740503e1ffc7047c562d7</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::EdgeData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rot</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ac0d833ca34bdab3452e9abf4c31822f7</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>div</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ac3f2039505ad380e35c592e44c3043ce</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, double alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0, int m=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>div</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a7dc894874cf6586a841a2dcfc748b64c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, double alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0, int m=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>div</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a58796e2e2d00d1968e265c89227029a7</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, double alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0, int m=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a323a2eab8129054a08424c0853ed7873</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, double alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a88033b813f3e6679d73904a31ab7efe7</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; dst, double alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a0cccc90b0521196a099b6b522111e1fa</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, double alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a9672ce98853eedcff9b9824a804865ec</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>grad</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ada3c304a2ad1dcff84acf58ee3376971</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a6d03768fecccdb355ea198661076256c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ac11b5155e46106d8159c06cc0e37291f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a5d0ae3a0556d6c58925d87a2061cb592</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>interp</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a88d36999905f40d7971ccb50ca9fbb1f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>laplace</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a80c18d9a36f8e08fa0b96aaf5fb63b8f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, double alpha, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, double gamma, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0, int m=0, int n=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>laplace</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a8bccaaf4a5815573f91ba2e89120eb6e</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, double alpha, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src1, double gamma, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0, int m=0, int n=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>laplace</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ad8b4b4529a5a667c082632b1aa38037f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; alpha, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, double gamma, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0, int m=0, int n=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>laplace</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ad602a2988d27105c19dba84888b458ab</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; alpha, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, double gamma, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0, int m=0, int n=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>vc_laplace</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a94639c26b65f5553e3467bd73bfb7373</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, double alpha, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; coef, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src1, double gamma, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int l=0, int m=0, int n=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a7a0c6b18ad01b1835e415aac73fe460a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, double alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>afaa6bdf034f047516975d9dde9c3d983</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0, int l=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a49bb3a25af5c04a76ff1162d7944d75f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0, int l=0, int m=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a0d90cda51faab692f175c26196fed37b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; dst, double alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>aa250cfd4cf2230429bcef3709d588bda</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0, int l=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a59bb75e34e272115a22a8b1f8a1d7a6d</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::FaceData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0, int l=0, int m=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>aa215f3a34574e6c4ef750b74f95ca1bc</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; dst, double alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ac134587f2fd3938b85cc255b49628524</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0, int l=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a56cb5e15fedfc60768ef56e1c2610f67</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; src1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0, int l=0, int m=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>aa15eff2d0462298293846377c1839be8</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, double alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>affdc64897c43742239cd902cdb1e4f01</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src1, double beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0, int l=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMultiply</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a95e358554939d402015510153d9465fa</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; alpha, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src1, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; beta, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; src2, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, int i=0, int j=0, int k=0, int l=0, int m=0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseL1Norm</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a7e87d204e00ec1d19d0d78ad6f7e2ef1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseL2Norm</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>aafa7b1e421edd6cc2b8ed2032f2aada5</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMaxNorm</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a61cd7193f65cbe4bc7686027ef37da89</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::CellData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseL1Norm</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>a22f7e1f55424a918cbde1d499eedb838</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseL2Norm</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>abe575b916cd62e2ac740d6cd4dedb330</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pointwiseMaxNorm</name>
      <anchorfile>class_i_b_t_k_1_1_patch_math_ops.html</anchorfile>
      <anchor>ab9abb6b9ac2e51b40cdb5f63c00224ff</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; dst, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::NodeData&lt; NDIM, double &gt; &gt; src, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScKrylovLinearSolver</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</filename>
    <base>IBTK::KrylovLinearSolver</base>
    <member kind="function">
      <type></type>
      <name>PETScKrylovLinearSolver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>ad238cee86d78cc8f2d99a69d25ded8dc</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix, MPI_Comm petsc_comm=PETSC_COMM_WORLD)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PETScKrylovLinearSolver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>a4d581fd57fbb226bb39a996ebc6fdb4b</anchor>
      <arglist>(const std::string &amp;object_name, const KSP &amp;petsc_ksp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PETScKrylovLinearSolver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>a5d7cb633ff53d315ed11ffa681d4d2ff</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setKSPType</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>a4d36f1d2e8f833f79b03250e198ea70a</anchor>
      <arglist>(const std::string &amp;ksp_type)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setOptionsPrefix</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>a469f57852e67fcb3711f858e1cd09ab8</anchor>
      <arglist>(const std::string &amp;options_prefix)</arglist>
    </member>
    <member kind="function">
      <type>const KSP &amp;</type>
      <name>getPETScKSP</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>a010c3112086e86b917ddc070350c46a0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; KrylovLinearSolver &gt;</type>
      <name>allocate_solver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>a0f8ffe44149f1147d37db7fbccb41860</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setOperator</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>a916b83d48d265389375c39716e384544</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LinearOperator &gt; A)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>aaad667b2a84c0fdd332cef6bf02b980d</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; LinearSolver &gt; pc_solver=NULL)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setNullspace</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>a8963f02086e155ff6ef96b0266ddb029</anchor>
      <arglist>(bool contains_constant_vec, const std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAIVectorReal_NDIM_double &gt; &gt; &amp;nullspace_basis_vecs=std::vector&lt; SAMRAI::tbox::Pointer&lt; SAMRAIVectorReal_NDIM_double &gt; &gt;())</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>a4a8a06cbd579a46a19d4ac8de14b4a0c</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>a43bf917ee0551640faab4a9f0b0c49b0</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_linear_solver.html</anchorfile>
      <anchor>a8fb9f6a1c4f1a15b87f7e961632be27d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScKrylovPoissonSolver</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_krylov_poisson_solver.html</filename>
    <base>IBTK::PETScKrylovLinearSolver</base>
    <base>IBTK::KrylovLinearSolverPoissonSolverInterface</base>
    <member kind="function">
      <type></type>
      <name>PETScKrylovPoissonSolver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_poisson_solver.html</anchorfile>
      <anchor>ac18b5df31a8a471058c11624055ab2db</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PETScKrylovPoissonSolver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_krylov_poisson_solver.html</anchorfile>
      <anchor>a10eaae779d023b5e55ac3ebb273a1398</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScLevelSolver</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</filename>
    <base>IBTK::LinearSolver</base>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>initializeSolverStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>abd46672f95b75680046ff35b413fe531</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)=0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>deallocateSolverStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a5d52aa0698168b0d32b062efb59e8183</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>copyToPETScVec</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a287a5008084baaf07eafec6d57ad2dec</anchor>
      <arglist>(Vec &amp;petsc_x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)=0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>copyFromPETScVec</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>af72ce52b6902d6886fa9945bb512e90a</anchor>
      <arglist>(Vec &amp;petsc_x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)=0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>setupKSPVecs</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a324df2fbb7b3f708a0173be231dbe044</anchor>
      <arglist>(Vec &amp;petsc_x, Vec &amp;petsc_b, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScMatLOWrapper</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_mat_l_o_wrapper.html</filename>
    <base>IBTK::LinearOperator</base>
    <member kind="function">
      <type></type>
      <name>PETScMatLOWrapper</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_l_o_wrapper.html</anchorfile>
      <anchor>a0a866960a50b72f9a6af7a140b0c46f5</anchor>
      <arglist>(const std::string &amp;object_name, const Mat &amp;petsc_mat)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PETScMatLOWrapper</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_l_o_wrapper.html</anchorfile>
      <anchor>aee0f31fcea61588a3f228a21ee591253</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Mat &amp;</type>
      <name>getPETScMat</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_l_o_wrapper.html</anchorfile>
      <anchor>a864e7f04e7873a0dc2e2a670bd2de603</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>apply</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_l_o_wrapper.html</anchorfile>
      <anchor>af071ff0d974c5c2e5ac891aeb0dc3259</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyAdd</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_l_o_wrapper.html</anchorfile>
      <anchor>a441dee7471961f605c20707f682e3e01</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;z)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_l_o_wrapper.html</anchorfile>
      <anchor>abf00529fb470bc4eb83686b99b522e1b</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_l_o_wrapper.html</anchorfile>
      <anchor>a967ae7f8a61adf5ef9b0fc7d360e085d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScMatUtilities</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_mat_utilities.html</filename>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>constructPatchLevelCCLaplaceOp</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_utilities.html</anchorfile>
      <anchor>a38233fff79b6e1d89c08a00ef15a15d9</anchor>
      <arglist>(Mat &amp;mat, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef, double data_time, const std::vector&lt; int &gt; &amp;num_dofs_per_proc, int dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>constructPatchLevelCCLaplaceOp</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_utilities.html</anchorfile>
      <anchor>a5b3eb1f81448e6440b0dfb8ac476dd2e</anchor>
      <arglist>(Mat &amp;mat, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, double data_time, const std::vector&lt; int &gt; &amp;num_dofs_per_proc, int dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>constructPatchLevelCCComplexLaplaceOp</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_utilities.html</anchorfile>
      <anchor>aaf99f063b3d87cea4f750d19decd57ae</anchor>
      <arglist>(Mat &amp;mat, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_real, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_imag, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef, double data_time, const std::vector&lt; int &gt; &amp;num_dofs_per_proc, int dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>constructPatchLevelCCComplexLaplaceOp</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_utilities.html</anchorfile>
      <anchor>ab930940a8eb7dbb85df60e1c19ee53fc</anchor>
      <arglist>(Mat &amp;mat, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_real, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_imag, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, double data_time, const std::vector&lt; int &gt; &amp;num_dofs_per_proc, int dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>constructPatchLevelSCLaplaceOp</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_utilities.html</anchorfile>
      <anchor>a3752503ca81938bb163d2d97bb9fc134</anchor>
      <arglist>(Mat &amp;mat, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, double data_time, const std::vector&lt; int &gt; &amp;num_dofs_per_proc, int dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>constructPatchLevelSCInterpOp</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_mat_utilities.html</anchorfile>
      <anchor>a5d880a98929ee6f9b832fe0033f475df</anchor>
      <arglist>(Mat &amp;mat, void(*interp_fcn)(double r_lower, double *w), int interp_stencil, Vec &amp;X_vec, const std::vector&lt; int &gt; &amp;num_dofs_per_proc, int dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScMFFDJacobianOperator</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_m_f_f_d_jacobian_operator.html</filename>
    <base>IBTK::JacobianOperator</base>
    <member kind="function">
      <type></type>
      <name>PETScMFFDJacobianOperator</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_m_f_f_d_jacobian_operator.html</anchorfile>
      <anchor>af24dfe53d4474614c9c8667866e79155</anchor>
      <arglist>(const std::string &amp;object_name, const std::string &amp;options_prefix=&quot;&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PETScMFFDJacobianOperator</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_m_f_f_d_jacobian_operator.html</anchorfile>
      <anchor>aed8379e4e7fbd9cc06eba0dda375d77e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setOperator</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_m_f_f_d_jacobian_operator.html</anchorfile>
      <anchor>a0fbf1d7a380d2d072603f354b7b2c24a</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; GeneralOperator &gt; F)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setNewtonKrylovSolver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_m_f_f_d_jacobian_operator.html</anchorfile>
      <anchor>ab345ee72e0a0e71138596e50ba531f04</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; PETScNewtonKrylovSolver &gt; nonlinear_solver)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>formJacobian</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_m_f_f_d_jacobian_operator.html</anchorfile>
      <anchor>a9180bb961c3d3c6e9982e4bd80b3e2a1</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;u)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt;</type>
      <name>getBaseVector</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_m_f_f_d_jacobian_operator.html</anchorfile>
      <anchor>a7f667afb3f3e60994831e083d7bdda51</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>apply</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_m_f_f_d_jacobian_operator.html</anchorfile>
      <anchor>ac784871244786e1be8cf8c97b822f808</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_m_f_f_d_jacobian_operator.html</anchorfile>
      <anchor>ad4332ba6108997970db79bb4a5678ada</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_m_f_f_d_jacobian_operator.html</anchorfile>
      <anchor>ac0b8b468945f386b64b300ab429cc117</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScNewtonKrylovSolver</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</filename>
    <base>IBTK::NewtonKrylovSolver</base>
    <member kind="function">
      <type></type>
      <name>PETScNewtonKrylovSolver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>a4c9888ec3a49b4a963f20a36285d93c4</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix, MPI_Comm petsc_comm=PETSC_COMM_WORLD)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PETScNewtonKrylovSolver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>a6b87ffe5bab80c996f6cf8819aa9da29</anchor>
      <arglist>(const std::string &amp;object_name, const SNES &amp;petsc_snes)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PETScNewtonKrylovSolver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>ad8db84cb8b7e3b2fe3107531bc9331ff</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setOptionsPrefix</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>a694674fa8214478f826d60ee7c3a62f9</anchor>
      <arglist>(const std::string &amp;options_prefix)</arglist>
    </member>
    <member kind="function">
      <type>const SNES &amp;</type>
      <name>getPETScSNES</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>a43332fbea2b2b0b3a8b77c1ad213dd87</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setOperator</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>a3e3afe3be1464fe15e23a8da9d32374b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; GeneralOperator &gt; op)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt;</type>
      <name>getSolutionVector</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>a111d7275bbb1af7e8d6e5da07dd810c7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt;</type>
      <name>getFunctionVector</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>ab7a1bee2711703ead94dba853d156a0e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setJacobian</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>a12172ea46196b2583cb3ddced038cc0f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; JacobianOperator &gt; J)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>ab51d1ffac3bcac2586b409c5ebf9a65c</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>a9bed200f6ff1e336ef3dfd4528ae148e</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>a0b719001075a5ffb946a517ef49006cb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; NewtonKrylovSolver &gt;</type>
      <name>allocate_solver</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_newton_krylov_solver.html</anchorfile>
      <anchor>a91dfff077d248ed653776f8f66d2f0b6</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScPCLSWrapper</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</filename>
    <base>IBTK::LinearSolver</base>
    <member kind="function">
      <type></type>
      <name>PETScPCLSWrapper</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a27d71163d6120133b6439398c665696a</anchor>
      <arglist>(const std::string &amp;object_name, const PC &amp;petsc_pc)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PETScPCLSWrapper</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a3108901e318698e01839d6f1068ac0f5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const PC &amp;</type>
      <name>getPETScPC</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a446fffbf6d63168710c7069757e6183b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a65a42a7f9fab2518284f5767a353fd4f</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a326e01517ac89232f5477bc160c6aa7f</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>aeda2460203abbac355502f45322eed50</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setInitialGuessNonzero</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>aeaec0254d26668ae210cffa2bfb2088b</anchor>
      <arglist>(bool initial_guess_nonzero=true)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>getInitialGuessNonzero</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a48a2472b694b6f3f5393c1466def315d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setMaxIterations</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a46871d043848c1d2c119f3841febaeaa</anchor>
      <arglist>(int max_iterations)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getMaxIterations</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a79a4454631cab17b3b4469dfb25dbb06</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setAbsoluteTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a009cedd8fa3ab8c12b6306f6d5590166</anchor>
      <arglist>(double abs_residual_tol)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getAbsoluteTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a18cb2bea908d886a77d22b4597d260f0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setRelativeTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a74e397a39b522f32dcd176ce59e1240d</anchor>
      <arglist>(double rel_residual_tol)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getRelativeTolerance</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a9b4a980a7d6c4bc212164a8ecd5a7a05</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumIterations</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>a20a3038f5744a0885da00d18794c5ca7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getResidualNorm</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_p_c_l_s_wrapper.html</anchorfile>
      <anchor>adae7d64284745dd3801d5988c031590b</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScSAMRAIVectorReal</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_s_a_m_r_a_i_vector_real.html</filename>
    <member kind="function" static="yes">
      <type>static Vec</type>
      <name>createPETScVector</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_a_m_r_a_i_vector_real.html</anchorfile>
      <anchor>ac740385b61425b896ff3741330dd3cb4</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, PetscScalar &gt; &gt; samrai_vec, MPI_Comm comm=PETSC_COMM_WORLD)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>destroyPETScVector</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_a_m_r_a_i_vector_real.html</anchorfile>
      <anchor>a9d0053102c380d3817ee94caab331708</anchor>
      <arglist>(Vec petsc_vec)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, PetscScalar &gt; &gt;</type>
      <name>getSAMRAIVector</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_a_m_r_a_i_vector_real.html</anchorfile>
      <anchor>a7141b19bff0002cf67b8b52817cc711d</anchor>
      <arglist>(Vec petsc_vec)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>replaceSAMRAIVector</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_a_m_r_a_i_vector_real.html</anchorfile>
      <anchor>a2e5164e1ac44fd9561554f9b49ea4254</anchor>
      <arglist>(Vec petsc_vec, SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, PetscScalar &gt; &gt; samrai_vec)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScSNESFunctionGOWrapper</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_function_g_o_wrapper.html</filename>
    <base>IBTK::GeneralOperator</base>
    <member kind="function">
      <type></type>
      <name>PETScSNESFunctionGOWrapper</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_function_g_o_wrapper.html</anchorfile>
      <anchor>a7b11a7ce57073e90d61d82f65afa5f24</anchor>
      <arglist>(const std::string &amp;object_name, const SNES &amp;petsc_snes, PetscErrorCode(*petsc_snes_form_func)(SNES, Vec, Vec, void *), void *petsc_snes_func_ctx)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PETScSNESFunctionGOWrapper</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_function_g_o_wrapper.html</anchorfile>
      <anchor>a55a889b110ef4086aa51e627c94d8771</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>apply</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_function_g_o_wrapper.html</anchorfile>
      <anchor>a67a01118d459da7a0cbfd39f84d74102</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_function_g_o_wrapper.html</anchorfile>
      <anchor>ab235536979dec1f12e0a7a8a780fa64e</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_function_g_o_wrapper.html</anchorfile>
      <anchor>ac3c30b1079cd19930c9ea8ba1ae212ad</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>PetscErrorCode(*)(SNES, Vec, Vec, void *)</type>
      <name>getPETScSNESFormFunction</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_function_g_o_wrapper.html</anchorfile>
      <anchor>a6a6c0ee9f2fc86ceb93988022428068b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const SNES &amp;</type>
      <name>getPETScSNES</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_function_g_o_wrapper.html</anchorfile>
      <anchor>ae838fa8963dfed5108d2833a3d04b1cc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>getPETScSNESFunctionContext</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_function_g_o_wrapper.html</anchorfile>
      <anchor>a9abcd0fceaf5b33ce658e88501c5afdb</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScSNESJacobianJOWrapper</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</filename>
    <base>IBTK::JacobianOperator</base>
    <member kind="function">
      <type></type>
      <name>PETScSNESJacobianJOWrapper</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</anchorfile>
      <anchor>ae5970284220265ea19d01f60fc0b9235</anchor>
      <arglist>(const std::string &amp;object_name, const SNES &amp;petsc_snes, PetscErrorCode(*petsc_snes_form_jac)(SNES, Vec, Mat *, Mat *, MatStructure *, void *), void *petsc_snes_jac_ctx)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PETScSNESJacobianJOWrapper</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</anchorfile>
      <anchor>a8aa1e6dbe62796de03153c65913909fe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>formJacobian</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</anchorfile>
      <anchor>aa55af546040b8989792261b2271a43cc</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &gt;</type>
      <name>getBaseVector</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</anchorfile>
      <anchor>a5dc792bacb7edf787d4eb6cb99a4df64</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>apply</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</anchorfile>
      <anchor>a0f0803df0c0191a90b17b3e4d06ceed5</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyAdd</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</anchorfile>
      <anchor>a50a2b11a6e48ef675b93ada3f2ee1521</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;z)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</anchorfile>
      <anchor>a523542faacb2fd02e834cedfed4b9cba</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</anchorfile>
      <anchor>aa7a144ebc49a1a5e38d81f18b2808b36</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>PetscErrorCode(*)(SNES, Vec, Mat *, Mat *, MatStructure *, void *)</type>
      <name>getPETScSNESFormJacobian</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</anchorfile>
      <anchor>afef15a582a7c58af26be4ead02a12c9c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const SNES &amp;</type>
      <name>getPETScSNES</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</anchorfile>
      <anchor>a1841cdbd63abeb91997a6d33cf470a90</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>getPETScSNESJacobianContext</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_s_n_e_s_jacobian_j_o_wrapper.html</anchorfile>
      <anchor>ac5ab99d5ab7f0f245128a75bbe045d31</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PETScVecUtilities</name>
    <filename>class_i_b_t_k_1_1_p_e_t_sc_vec_utilities.html</filename>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>copyToPatchLevelVec</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_vec_utilities.html</anchorfile>
      <anchor>a20630e823e550e40ab9d4409a5238f9e</anchor>
      <arglist>(Vec &amp;vec, int data_idx, int dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>copyFromPatchLevelVec</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_vec_utilities.html</anchorfile>
      <anchor>a635d40073d1cea7c22d44c0f066e2d7a</anchor>
      <arglist>(Vec &amp;vec, int data_idx, int dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; data_synch_sched, SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt; ghost_fill_sched)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt;</type>
      <name>constructDataSynchSchedule</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_vec_utilities.html</anchorfile>
      <anchor>a2a7e87ec3237c2957af8f31ab2254ce2</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; SAMRAI::xfer::RefineSchedule&lt; NDIM &gt; &gt;</type>
      <name>constructGhostFillSchedule</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_vec_utilities.html</anchorfile>
      <anchor>a9c8f86ab80696952cf1122a7f55cb999</anchor>
      <arglist>(int data_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>constructPatchLevelDOFIndices</name>
      <anchorfile>class_i_b_t_k_1_1_p_e_t_sc_vec_utilities.html</anchorfile>
      <anchor>aec166bd294a3a133d5c24bfc015f9c62</anchor>
      <arglist>(std::vector&lt; int &gt; &amp;num_dofs_per_proc, int dof_index_idx, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PhysicalBoundaryUtilities</name>
    <filename>class_i_b_t_k_1_1_physical_boundary_utilities.html</filename>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>isLower</name>
      <anchorfile>class_i_b_t_k_1_1_physical_boundary_utilities.html</anchorfile>
      <anchor>a58c9e6bd0bbdf33a7889e67d985ed830</anchor>
      <arglist>(int loc, int codim, int direction)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>isUpper</name>
      <anchorfile>class_i_b_t_k_1_1_physical_boundary_utilities.html</anchorfile>
      <anchor>ac6a0132a1e5ccd1be5671650f76b8361</anchor>
      <arglist>(int loc, int codim, int direction)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Array&lt; SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &gt;</type>
      <name>getPhysicalBoundaryCodim1Boxes</name>
      <anchorfile>class_i_b_t_k_1_1_physical_boundary_utilities.html</anchorfile>
      <anchor>a46b91760436cb92dd143ae6d70af27d1</anchor>
      <arglist>(const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Array&lt; SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &gt;</type>
      <name>getPhysicalBoundaryCodim2Boxes</name>
      <anchorfile>class_i_b_t_k_1_1_physical_boundary_utilities.html</anchorfile>
      <anchor>a50fb2dc7cbf412360999125650707910</anchor>
      <arglist>(const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Array&lt; SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &gt;</type>
      <name>getPhysicalBoundaryCodim3Boxes</name>
      <anchorfile>class_i_b_t_k_1_1_physical_boundary_utilities.html</anchorfile>
      <anchor>a6c3bd0f6bbfe726f97e1df0a2275ec8c</anchor>
      <arglist>(const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::hier::BoundaryBox&lt; NDIM &gt;</type>
      <name>trimBoundaryCodim1Box</name>
      <anchorfile>class_i_b_t_k_1_1_physical_boundary_utilities.html</anchorfile>
      <anchor>a06a8e9bd0b70df1be4699fc118909580</anchor>
      <arglist>(const SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &amp;bdry_box, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::hier::Box&lt; NDIM &gt;</type>
      <name>makeSideBoundaryCodim1Box</name>
      <anchorfile>class_i_b_t_k_1_1_physical_boundary_utilities.html</anchorfile>
      <anchor>a125a3b7858ded6d956ea567598dd0471</anchor>
      <arglist>(const SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &amp;boundary_box)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PoissonFACPreconditioner</name>
    <filename>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner.html</filename>
    <base>IBTK::FACPreconditioner</base>
    <base>IBTK::PoissonSolver</base>
    <member kind="function">
      <type></type>
      <name>PoissonFACPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner.html</anchorfile>
      <anchor>a272a73d6149ee5cfb64e87d6bbf9d96c</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; PoissonFACPreconditionerStrategy &gt; fac_strategy, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PoissonFACPreconditioner</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner.html</anchorfile>
      <anchor>ad63a8a9e615467066c044c5fc6578910</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPoissonSpecifications</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner.html</anchorfile>
      <anchor>ab70396c7a6e3735cf12c520365f142ca</anchor>
      <arglist>(const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoef</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner.html</anchorfile>
      <anchor>a0f01a5dedbc90cad76c495dc865ced5c</anchor>
      <arglist>(SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBcCoefs</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner.html</anchorfile>
      <anchor>ab9f80860194d3992bd3371d62eba4c2f</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PoissonFACPreconditionerStrategy</name>
    <filename>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</filename>
    <base>IBTK::FACPreconditionerStrategy</base>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setSmootherType</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a1f3580140a477cac1dab8e0dfa6e3519</anchor>
      <arglist>(const std::string &amp;smoother_type)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setCoarseSolverType</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>ab7840ccbce1ea22167193faeee4b56cd</anchor>
      <arglist>(const std::string &amp;coarse_solver_type)=0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>initializeOperatorStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>aadc2c0ce6a0f2877c24e2d802f68c750</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs, int coarsest_reset_ln, int finest_reset_ln)=0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>deallocateOperatorStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_f_a_c_preconditioner_strategy.html</anchorfile>
      <anchor>a613fdc4b79eb2d7bcbe44849495f9b0a</anchor>
      <arglist>(int coarsest_reset_ln, int finest_reset_ln)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::PoissonSolver</name>
    <filename>class_i_b_t_k_1_1_poisson_solver.html</filename>
    <base virtualness="virtual">IBTK::GeneralSolver</base>
  </compound>
  <compound kind="class">
    <name>IBTK::PoissonUtilities</name>
    <filename>class_i_b_t_k_1_1_poisson_utilities.html</filename>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>computeCCMatrixCoefficients</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_utilities.html</anchorfile>
      <anchor>ac3a4b9fe2a4d23e4761555b6e5232610</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;matrix_coefficients, const std::vector&lt; SAMRAI::hier::Index&lt; NDIM &gt; &gt; &amp;stencil, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef, double data_time)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>computeCCMatrixCoefficients</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_utilities.html</anchorfile>
      <anchor>a9e89b05cc43b0e53940dde219be57eae</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;matrix_coefficients, const std::vector&lt; SAMRAI::hier::Index&lt; NDIM &gt; &gt; &amp;stencil, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, double data_time)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>computeCCComplexMatrixCoefficients</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_utilities.html</anchorfile>
      <anchor>a8c0a9a35a9ba24ae674d69bcfe49ea4c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;matrix_coefficients, const std::vector&lt; SAMRAI::hier::Index&lt; NDIM &gt; &gt; &amp;stencil, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_real, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_imag, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef, double data_time)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>computeCCComplexMatrixCoefficients</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_utilities.html</anchorfile>
      <anchor>ad6908b00ebe09d6c92ca3c51d062ffbd</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;matrix_coefficients, const std::vector&lt; SAMRAI::hier::Index&lt; NDIM &gt; &gt; &amp;stencil, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_real, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_imag, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, double data_time)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>computeSCMatrixCoefficients</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_utilities.html</anchorfile>
      <anchor>aed6148e0df3f4ffc00c617eed7f5fe6f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, SAMRAI::pdat::SideData&lt; NDIM, double &gt; &amp;matrix_coefficients, const std::vector&lt; SAMRAI::hier::Index&lt; NDIM &gt; &gt; &amp;stencil, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, double data_time)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>adjustCCBoundaryRhsEntries</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_utilities.html</anchorfile>
      <anchor>a641b42fd21806546ab904d9939ed3e18</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;rhs_data, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef, double data_time, bool homogeneous_bc)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>adjustCCBoundaryRhsEntries</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_utilities.html</anchorfile>
      <anchor>a8ee64d7bcd26fdac1f12221f12c46cb5</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;rhs_data, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, double data_time, bool homogeneous_bc)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>adjustCCComplexBoundaryRhsEntries</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_utilities.html</anchorfile>
      <anchor>afb7a7391ea526431834afbb0482cbc69</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;rhs_data, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_real, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_imag, SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; *bc_coef, double data_time, bool homogeneous_bc)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>adjustCCComplexBoundaryRhsEntries</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_utilities.html</anchorfile>
      <anchor>a4ed5a70e8685b0b6ed2e0c9b213f4601</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, SAMRAI::pdat::CellData&lt; NDIM, double &gt; &amp;rhs_data, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_real, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec_imag, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, double data_time, bool homogeneous_bc)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>adjustSCBoundaryRhsEntries</name>
      <anchorfile>class_i_b_t_k_1_1_poisson_utilities.html</anchorfile>
      <anchor>ab6cfd5f4c19cec4f70540051d7c4e5f1</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, SAMRAI::pdat::SideData&lt; NDIM, double &gt; &amp;rhs_data, const SAMRAI::solv::PoissonSpecifications &amp;poisson_spec, const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;bc_coefs, double data_time, bool homogeneous_bc)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::RefinePatchStrategySet</name>
    <filename>class_i_b_t_k_1_1_refine_patch_strategy_set.html</filename>
    <base>RefinePatchStrategy&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>RefinePatchStrategySet</name>
      <anchorfile>class_i_b_t_k_1_1_refine_patch_strategy_set.html</anchorfile>
      <anchor>a041a1f7b1f8f6f24f85c4ee2ecc5cb58</anchor>
      <arglist>(InputIterator first, InputIterator last, bool managed=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RefinePatchStrategySet</name>
      <anchorfile>class_i_b_t_k_1_1_refine_patch_strategy_set.html</anchorfile>
      <anchor>af5117af7875023e98e074186a47e663f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchorfile>class_i_b_t_k_1_1_refine_patch_strategy_set.html</anchorfile>
      <anchor>a5c216768f51c0bf4552ee5fbc7b2b80a</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt;</type>
      <name>getRefineOpStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_refine_patch_strategy_set.html</anchorfile>
      <anchor>accc287047b6fb2b4d0784c7dbd60ee37</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_refine_patch_strategy_set.html</anchorfile>
      <anchor>a799c365dbde692d7b9901ebfa5f94ab5</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessRefine</name>
      <anchorfile>class_i_b_t_k_1_1_refine_patch_strategy_set.html</anchorfile>
      <anchor>a5dbdaa94fb9e1511462754831d4d2c40</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;fine_box, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>preprocessRefineBoxes</name>
      <anchorfile>class_i_b_t_k_1_1_refine_patch_strategy_set.html</anchorfile>
      <anchor>a984c7e0a37e957121e5de454ba5ac337</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::BoxList&lt; NDIM &gt; &amp;fine_boxes, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>postprocessRefineBoxes</name>
      <anchorfile>class_i_b_t_k_1_1_refine_patch_strategy_set.html</anchorfile>
      <anchor>a012b30fd9ee2b665f948064dccb9bc8c</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;fine, const SAMRAI::hier::Patch&lt; NDIM &gt; &amp;coarse, const SAMRAI::hier::BoxList&lt; NDIM &gt; &amp;fine_boxes, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ratio)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::RobinPhysBdryPatchStrategy</name>
    <filename>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</filename>
    <base>RefinePatchStrategy&lt; NDIM &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>accumulateFromPhysicalBoundaryData</name>
      <anchorfile>class_i_b_t_k_1_1_robin_phys_bdry_patch_strategy.html</anchorfile>
      <anchor>ad7c9ec034342035573d709a0c271455a</anchor>
      <arglist>(SAMRAI::hier::Patch&lt; NDIM &gt; &amp;patch, double fill_time, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::SCLaplaceOperator</name>
    <filename>class_i_b_t_k_1_1_s_c_laplace_operator.html</filename>
    <base>IBTK::LaplaceOperator</base>
    <member kind="function">
      <type></type>
      <name>SCLaplaceOperator</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_laplace_operator.html</anchorfile>
      <anchor>a07bd0c8d177e88b2b2df072af84de035</anchor>
      <arglist>(const std::string &amp;object_name, bool homogeneous_bc=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SCLaplaceOperator</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_laplace_operator.html</anchorfile>
      <anchor>ae859c4fd0feed7788808352c7c5967e1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>apply</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_laplace_operator.html</anchorfile>
      <anchor>ab9cedd783e82b80168d11470a29f7f10</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_laplace_operator.html</anchorfile>
      <anchor>a8d32036c5e4fa5235d700e65cfdde618</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;in, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_laplace_operator.html</anchorfile>
      <anchor>a4e8904553c2f48f0887923a35b77128a</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::SCPoissonHypreLevelSolver</name>
    <filename>class_i_b_t_k_1_1_s_c_poisson_hypre_level_solver.html</filename>
    <base>IBTK::LinearSolver</base>
    <base>IBTK::PoissonSolver</base>
    <member kind="function">
      <type></type>
      <name>SCPoissonHypreLevelSolver</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>abfabf414f81e736e3be9c66dc77b3249</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SCPoissonHypreLevelSolver</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>a120b4d80a2e39618f0d8efb57fcfd093</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveSystem</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>ad8895163491905cfcc8a39e841f06370</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>a3d5bd3982d29cc498b5a0c170b2b93c8</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateSolverState</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>af09e127222639a34fbbd43a609bccfc7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;</type>
      <name>allocate_solver</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_hypre_level_solver.html</anchorfile>
      <anchor>ae4364590a8165ed6b69ad8e5a88dd93d</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::SCPoissonPETScLevelSolver</name>
    <filename>class_i_b_t_k_1_1_s_c_poisson_p_e_t_sc_level_solver.html</filename>
    <base>IBTK::PETScLevelSolver</base>
    <base>IBTK::PoissonSolver</base>
    <member kind="function">
      <type></type>
      <name>SCPoissonPETScLevelSolver</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a004bf9f0f98307916b3b28095c445880</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SCPoissonPETScLevelSolver</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a781a53bf52b9629f52057bc746e2bb83</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;</type>
      <name>allocate_solver</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>ac4beae7249cac8c96661a702af6f39ec</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initializeSolverStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>ae38d56f60a3b0da2b83a75b036a3e2e7</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>deallocateSolverStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a930aa1ff01ef78dc12b0cc39aa23dbec</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>copyToPETScVec</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a9b07a0dd8d53c51d45763c744bb33f04</anchor>
      <arglist>(Vec &amp;petsc_x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>copyFromPETScVec</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a2d1ca903772e453da956644c2df07b55</anchor>
      <arglist>(Vec &amp;petsc_x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>setupKSPVecs</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_p_e_t_sc_level_solver.html</anchorfile>
      <anchor>a3906e7c52275dcbcaa5ad5575ba1e5a9</anchor>
      <arglist>(Vec &amp;petsc_x, Vec &amp;petsc_b, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;x, SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;b, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; patch_level)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::SCPoissonPointRelaxationFACOperator</name>
    <filename>class_i_b_t_k_1_1_s_c_poisson_point_relaxation_f_a_c_operator.html</filename>
    <base>IBTK::PoissonFACPreconditionerStrategy</base>
    <member kind="function">
      <type></type>
      <name>SCPoissonPointRelaxationFACOperator</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>aa55af0b7555b5bebb897c7df88f724de</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SCPoissonPointRelaxationFACOperator</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a27dbfa75b9aa318655d25282de56d07f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setSmootherType</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>aec88da4111e68d1c4be9e0656842e894</anchor>
      <arglist>(const std::string &amp;smoother_type)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setCoarseSolverType</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a53e370f99d3a3d52766b6bd1de875d3b</anchor>
      <arglist>(const std::string &amp;coarse_solver_type)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>smoothError</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>ae4b5a37f93de5593e08a3320dc811fd1</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;error, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, int level_num, int num_sweeps, bool performing_pre_sweeps, bool performing_post_sweeps)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>solveCoarsestLevel</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a32a036b94a9f82f2f6096248d99cf0b3</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;error, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, int coarsest_ln)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>computeResidual</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>a3f7126771751d3d5877683533572b038</anchor>
      <arglist>(SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;residual, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs, int coarsest_level_num, int finest_level_num)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;</type>
      <name>allocate_solver</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>ad53b186828197c593b8575015cf5f18a</anchor>
      <arglist>(const std::string &amp;object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; input_db, const std::string &amp;default_options_prefix)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initializeOperatorStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>ae0dc675aa419b49217620d49cd1c27df</anchor>
      <arglist>(const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;solution, const SAMRAI::solv::SAMRAIVectorReal&lt; NDIM, double &gt; &amp;rhs, int coarsest_reset_ln, int finest_reset_ln)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>deallocateOperatorStateSpecialized</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_point_relaxation_f_a_c_operator.html</anchorfile>
      <anchor>ab7ea197822bab74c4225ad8ae5602bdf</anchor>
      <arglist>(int coarsest_reset_ln, int finest_reset_ln)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::SCPoissonSolverManager</name>
    <filename>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</filename>
    <member kind="typedef">
      <type>SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;(*</type>
      <name>SolverMaker</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>af7811ffb97edd7d759eda87021e57d6f</anchor>
      <arglist>)(const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;</type>
      <name>allocateSolver</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>a58b62d166b5b7f06f0d373603bf44391</anchor>
      <arglist>(const std::string &amp;solver_type, const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; PoissonSolver &gt;</type>
      <name>allocateSolver</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>a82ce7a8500eab7540cc0c28529bcd6ff</anchor>
      <arglist>(const std::string &amp;solver_type, const std::string &amp;solver_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; solver_input_db, const std::string &amp;solver_default_options_prefix, const std::string &amp;precond_type, const std::string &amp;precond_object_name, SAMRAI::tbox::Pointer&lt; SAMRAI::tbox::Database &gt; precond_input_db, const std::string &amp;precond_default_options_prefix) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>registerSolverFactoryFunction</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>abd58d60f93bd06789811b9df7c94c830</anchor>
      <arglist>(const std::string &amp;solver_type, SolverMaker solver_maker)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static SCPoissonSolverManager *</type>
      <name>getManager</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>a77a14da09060d12626cd2b672a15d17b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>freeManager</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>af17df459bf6dbc90925ee2babba0db22</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>UNDEFINED</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>a78a581f5cdc498c0c3f005ad67076c06</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT_KRYLOV_SOLVER</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>a734f4e290584e7921cc06189b3d7b729</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT_FAC_PRECONDITIONER</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>a5d52fef4b81a140e66998fdec933a0f6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const std::string</type>
      <name>DEFAULT_LEVEL_SOLVER</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>a926278877661cd22150863c56ff9ee42</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>SCPoissonSolverManager</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>ab1bf6811ee18c47bafad680785cbacbd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>~SCPoissonSolverManager</name>
      <anchorfile>class_i_b_t_k_1_1_s_c_poisson_solver_manager.html</anchorfile>
      <anchor>a753cf31e5929ec4801d26096189babdb</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::SideDataSynchronization</name>
    <filename>class_i_b_t_k_1_1_side_data_synchronization.html</filename>
    <class kind="class">IBTK::SideDataSynchronization::SynchronizationTransactionComponent</class>
    <member kind="function">
      <type></type>
      <name>SideDataSynchronization</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization.html</anchorfile>
      <anchor>a565a43ea387ba76a0b096d4b91dd7b2f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SideDataSynchronization</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization.html</anchorfile>
      <anchor>ab5c0e2d80e4fbe4035193488add8739c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization.html</anchorfile>
      <anchor>a110e3a1fdd889cdeb757b2ad52b519b6</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;transaction_comp, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization.html</anchorfile>
      <anchor>a87479707d9906ff0912f7cdd12b7619a</anchor>
      <arglist>(const std::vector&lt; SynchronizationTransactionComponent &gt; &amp;transaction_comps, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization.html</anchorfile>
      <anchor>acd85d1993a2ec78232a1f9ee9b0fc47b</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;transaction_comps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTransactionComponents</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization.html</anchorfile>
      <anchor>a3416d2dab105a28bfe7f9d05e7761e22</anchor>
      <arglist>(const std::vector&lt; SynchronizationTransactionComponent &gt; &amp;transaction_comps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deallocateOperatorState</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization.html</anchorfile>
      <anchor>a26b37095f461121aa3751edd7c2dad1e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>synchronizeData</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization.html</anchorfile>
      <anchor>ac18badd68d1e83a60c708cd92d4d380e</anchor>
      <arglist>(double fill_time)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::SideDataSynchronization::SynchronizationTransactionComponent</name>
    <filename>class_i_b_t_k_1_1_side_data_synchronization_1_1_synchronization_transaction_component.html</filename>
    <member kind="function">
      <type></type>
      <name>SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>a596206277010a570c7a8e19c5d0c1272</anchor>
      <arglist>(int data_idx=-1, const std::string &amp;coarsen_op_name=&quot;NONE&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>a4cc523c68c5b5cb12cc08897869ec44a</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;from)</arglist>
    </member>
    <member kind="function">
      <type>SynchronizationTransactionComponent &amp;</type>
      <name>operator=</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>ab82ba2e919a613f11653fdfbc9c08a37</anchor>
      <arglist>(const SynchronizationTransactionComponent &amp;that)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SynchronizationTransactionComponent</name>
      <anchorfile>class_i_b_t_k_1_1_side_data_synchronization_1_1_synchronization_transaction_component.html</anchorfile>
      <anchor>aa3129c0231935fded21bed422a13c40c</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::SideNoCornersFillPattern</name>
    <filename>class_i_b_t_k_1_1_side_no_corners_fill_pattern.html</filename>
    <base>VariableFillPattern&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>SideNoCornersFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_side_no_corners_fill_pattern.html</anchorfile>
      <anchor>a40c79a8fb75c8abf08b48496abbe8186</anchor>
      <arglist>(int stencil_width, bool include_dst_patch_box, bool include_edges_on_dst_level, bool include_edges_on_src_level)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SideNoCornersFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_side_no_corners_fill_pattern.html</anchorfile>
      <anchor>ad3383fe487a71b4e82ab4dfcdeef7307</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BoxOverlap&lt; NDIM &gt; &gt;</type>
      <name>calculateOverlap</name>
      <anchorfile>class_i_b_t_k_1_1_side_no_corners_fill_pattern.html</anchorfile>
      <anchor>afba912620ac42f066ddeae51621ef249</anchor>
      <arglist>(const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;dst_geometry, const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;src_geometry, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;dst_patch_box, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;src_mask, bool overwrite_interior, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BoxOverlap&lt; NDIM &gt; &gt;</type>
      <name>calculateOverlapOnLevel</name>
      <anchorfile>class_i_b_t_k_1_1_side_no_corners_fill_pattern.html</anchorfile>
      <anchor>a8bf69c8cbd8aa6188e8b5a085d0ec8ee</anchor>
      <arglist>(const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;dst_geometry, const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;src_geometry, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;dst_patch_box, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;src_mask, bool overwrite_interior, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset, int dst_level_num, int src_level_num) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setTargetPatchLevelNumber</name>
      <anchorfile>class_i_b_t_k_1_1_side_no_corners_fill_pattern.html</anchorfile>
      <anchor>aa4bf2ea48bcf4237ecc818c710a2c589</anchor>
      <arglist>(int level_num)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_side_no_corners_fill_pattern.html</anchorfile>
      <anchor>a4d36d409720213cb2874d8626177336a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getPatternName</name>
      <anchorfile>class_i_b_t_k_1_1_side_no_corners_fill_pattern.html</anchorfile>
      <anchor>aa6aca7db28747653ec2d8f7b40b0298b</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::SideSynchCopyFillPattern</name>
    <filename>class_i_b_t_k_1_1_side_synch_copy_fill_pattern.html</filename>
    <base>VariableFillPattern&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>SideSynchCopyFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_side_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a758831abec1f0d62e770e3033be739ba</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SideSynchCopyFillPattern</name>
      <anchorfile>class_i_b_t_k_1_1_side_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a8c5ee5d36528d9c28b8a6ff3706c4d26</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BoxOverlap&lt; NDIM &gt; &gt;</type>
      <name>calculateOverlap</name>
      <anchorfile>class_i_b_t_k_1_1_side_synch_copy_fill_pattern.html</anchorfile>
      <anchor>abc19d09e40f7d233ccfbb6753b84ea2c</anchor>
      <arglist>(const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;dst_geometry, const SAMRAI::hier::BoxGeometry&lt; NDIM &gt; &amp;src_geometry, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;dst_patch_box, const SAMRAI::hier::Box&lt; NDIM &gt; &amp;src_mask, bool overwrite_interior, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;src_offset) const </arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;</type>
      <name>getStencilWidth</name>
      <anchorfile>class_i_b_t_k_1_1_side_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a5d950d8c1be94ab64161abc4e212f24e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::string &amp;</type>
      <name>getPatternName</name>
      <anchorfile>class_i_b_t_k_1_1_side_synch_copy_fill_pattern.html</anchorfile>
      <anchor>a693a34aaad66bfc67117fefd0b5be44e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::StaggeredPhysicalBoundaryHelper</name>
    <filename>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</filename>
    <member kind="function">
      <type></type>
      <name>StaggeredPhysicalBoundaryHelper</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>a714856098ca67a3320435fe727add4d6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StaggeredPhysicalBoundaryHelper</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>aded6ae60b8d1c1929dcb2a2b4e21fd29</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>copyDataAtDirichletBoundaries</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>a22a9ff6a102fdc31707c5e8ac0a0e3ab</anchor>
      <arglist>(int u_out_data_idx, int u_in_data_idx, int coarsest_ln=-1, int finest_ln=-1) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>copyDataAtDirichletBoundaries</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>a686a4f2b950fda9a4e80a2967744ee09</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; u_out_data, SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, double &gt; &gt; u_in_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setupMaskingFunction</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>a4d2da84edaed424103db2e275b2ac596</anchor>
      <arglist>(int mask_data_idx, int coarsest_ln=-1, int finest_ln=-1) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setupMaskingFunction</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>a6b793dc1cfedb1609e7b05e296a1c6de</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::pdat::SideData&lt; NDIM, int &gt; &gt; u_data, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>patchTouchesDirichletBoundary</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>a623a169aab8ef0e9f38dd80dd372d19f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>patchTouchesDirichletBoundaryAxis</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>a6e8db9d57b97ee2410b0bde552010562</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch, const unsigned int axis) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cacheBcCoefData</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>a5d368884958c2d72a8a5a52840a879b7</anchor>
      <arglist>(const std::vector&lt; SAMRAI::solv::RobinBcCoefStrategy&lt; NDIM &gt; * &gt; &amp;u_bc_coefs, double fill_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>clearBcCoefData</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>aecb92d79d0d734bbb25c634c6273c8ed</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected" static="yes">
      <type>static void</type>
      <name>setupBcCoefBoxes</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>a1827cc2e92cf3282cd81c290e31cd873</anchor>
      <arglist>(SAMRAI::hier::Box&lt; NDIM &gt; &amp;bc_coef_box, SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &amp;trimmed_bdry_box, const SAMRAI::hier::BoundaryBox&lt; NDIM &gt; &amp;bdry_box, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::Patch&lt; NDIM &gt; &gt; patch)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt;</type>
      <name>d_hierarchy</name>
      <anchorfile>class_i_b_t_k_1_1_staggered_physical_boundary_helper.html</anchorfile>
      <anchor>a47cd4028cd0ac69a36c745e3f43866d4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::StandardTagAndInitStrategySet</name>
    <filename>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</filename>
    <base>StandardTagAndInitStrategy&lt; NDIM &gt;</base>
    <member kind="function">
      <type></type>
      <name>StandardTagAndInitStrategySet</name>
      <anchorfile>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</anchorfile>
      <anchor>a0bd1d748ceb96053f8ca306c3cdcca13</anchor>
      <arglist>(InputIterator first, InputIterator last, bool managed=true)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StandardTagAndInitStrategySet</name>
      <anchorfile>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</anchorfile>
      <anchor>a4cf967e3f53dc9347af502852f14434c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>getLevelDt</name>
      <anchorfile>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</anchorfile>
      <anchor>a24722ec4e2e58ad2d36c6c28e69ea7da</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; level, double dt_time, bool initial_time)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>advanceLevel</name>
      <anchorfile>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</anchorfile>
      <anchor>a97274e664600e3fa05f2550ea115bd49</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; level, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, double current_time, double new_time, bool first_step, bool last_step, bool regrid_advance=false)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetTimeDependentData</name>
      <anchorfile>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</anchorfile>
      <anchor>adebe4da5cb7a696d0d331c7e51629a2b</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; level, double new_time, bool can_be_refined)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetDataToPreadvanceState</name>
      <anchorfile>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</anchorfile>
      <anchor>aa74c9b126117e7aa58c9cdd8c0e77ca8</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; level)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeLevelData</name>
      <anchorfile>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</anchorfile>
      <anchor>a2e8b5bbd917f88b5ce6a79674ac8aa13</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double init_data_time, bool can_be_refined, bool initial_time, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt; old_level=SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchLevel&lt; NDIM &gt; &gt;(NULL), bool allocate_data=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetHierarchyConfiguration</name>
      <anchorfile>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</anchorfile>
      <anchor>a4bf3834cbcecde8caa4921a7749c1013</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyGradientDetector</name>
      <anchorfile>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</anchorfile>
      <anchor>a2eeae6cfd174cdcb06aac92c65e3f01f</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::BasePatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, double error_data_time, int tag_index, bool initial_time, bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyRichardsonExtrapolation</name>
      <anchorfile>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</anchorfile>
      <anchor>adfa048eb5256fd74b341f4f516ea9483</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; level, double error_data_time, int tag_index, double deltat, int error_coarsen_ratio, bool initial_time, bool uses_gradient_detector_too)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>coarsenDataForRichardsonExtrapolation</name>
      <anchorfile>class_i_b_t_k_1_1_standard_tag_and_init_strategy_set.html</anchorfile>
      <anchor>a8bf2665d426ba2aeeb8004343359282c</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchHierarchy&lt; NDIM &gt; &gt; hierarchy, int level_number, SAMRAI::tbox::Pointer&lt; SAMRAI::hier::PatchLevel&lt; NDIM &gt; &gt; coarser_level, double coarsen_data_time, bool before_advance)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::Streamable</name>
    <filename>class_i_b_t_k_1_1_streamable.html</filename>
    <member kind="function">
      <type></type>
      <name>Streamable</name>
      <anchorfile>class_i_b_t_k_1_1_streamable.html</anchorfile>
      <anchor>ad622395dbcdf1d48d63bd5a6c14094c4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Streamable</name>
      <anchorfile>class_i_b_t_k_1_1_streamable.html</anchorfile>
      <anchor>a5572cdd9379c74e3c82218cbb66afd14</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getStreamableClassID</name>
      <anchorfile>class_i_b_t_k_1_1_streamable.html</anchorfile>
      <anchor>a53009b530b448a9b5ae8a63d03cf7e67</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_t_k_1_1_streamable.html</anchorfile>
      <anchor>a283af2eec5511a777e86ff47dfa5c690</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_t_k_1_1_streamable.html</anchorfile>
      <anchor>ac2ed603d8ebcb2b158e5b3e4a18ec61b</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerPeriodicShift</name>
      <anchorfile>class_i_b_t_k_1_1_streamable.html</anchorfile>
      <anchor>adc3aac0c985d54e0ca23014625fb8a15</anchor>
      <arglist>(const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset, const Vector &amp;displacement)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::StreamableFactory</name>
    <filename>class_i_b_t_k_1_1_streamable_factory.html</filename>
    <member kind="function">
      <type></type>
      <name>StreamableFactory</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_factory.html</anchorfile>
      <anchor>a49a6a57adf2b32832674507629ad667c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~StreamableFactory</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_factory.html</anchorfile>
      <anchor>a735c25d13917ca5c6fb681b227f48364</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getStreamableClassID</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_factory.html</anchorfile>
      <anchor>a854bc32e9c1e8f216cbd914eb12b2fdc</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setStreamableClassID</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_factory.html</anchorfile>
      <anchor>a24265196f6417863627846ed39721582</anchor>
      <arglist>(int class_id)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual SAMRAI::tbox::Pointer&lt; Streamable &gt;</type>
      <name>unpackStream</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_factory.html</anchorfile>
      <anchor>a521be10626d985f3a9d37694647f27ab</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IBTK::StreamableManager</name>
    <filename>class_i_b_t_k_1_1_streamable_manager.html</filename>
    <member kind="function">
      <type>bool</type>
      <name>checkFactoryRegistration</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>ae877253386b6ea195c12e453a909c9bd</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; StreamableFactory &gt; factory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>registerFactory</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>a52b81b17cbec98b5a7e0239debf775e4</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; StreamableFactory &gt; factory)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>a416e0e284ff922e8864a394d4bc10191</anchor>
      <arglist>(SAMRAI::tbox::Pointer&lt; Streamable &gt; data_item) const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getDataStreamSize</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>a2f4ede28d7c6e6aacbf4128de0561882</anchor>
      <arglist>(const std::vector&lt; SAMRAI::tbox::Pointer&lt; Streamable &gt; &gt; &amp;data_items) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>ac1d42dff9916a756d15f18a91463128e</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, SAMRAI::tbox::Pointer&lt; Streamable &gt; data_item)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>packStream</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>a89f9fd4bea4243be356edb548788355b</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, std::vector&lt; SAMRAI::tbox::Pointer&lt; Streamable &gt; &gt; &amp;data_items)</arglist>
    </member>
    <member kind="function">
      <type>SAMRAI::tbox::Pointer&lt; Streamable &gt;</type>
      <name>unpackStream</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>a1fc50986146be1bfaee73a88f4384cf5</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unpackStream</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>a8e73535f3e7ee635e61a5b52200ecada</anchor>
      <arglist>(SAMRAI::tbox::AbstractStream &amp;stream, const SAMRAI::hier::IntVector&lt; NDIM &gt; &amp;offset, std::vector&lt; SAMRAI::tbox::Pointer&lt; Streamable &gt; &gt; &amp;data_items)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static StreamableManager *</type>
      <name>getManager</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>a52c10d437abab763883379566c632c3b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>freeManager</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>ad6a8c641d0643c17c81f862f511da98e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>getUnregisteredID</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>a857ad316d177a284ec663cc80f56edb3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>StreamableManager</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>a955bb029ae3c73bd3b8802084ad771e5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>~StreamableManager</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>a08d0faa65286c77aee9e9924a2c68349</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected" static="yes">
      <type>static int</type>
      <name>createUniqueID</name>
      <anchorfile>class_i_b_t_k_1_1_streamable_manager.html</anchorfile>
      <anchor>ab464e5ed1646efccc27078fa16ac3a2a</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/boundary</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/boundary/</path>
    <filename>dir_23d86e090902a6183d7a47e1c008559c.html</filename>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/boundary/cf_interface</dir>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/boundary/physical_boundary</dir>
    <file>HierarchyGhostCellInterpolation.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/boundary/cf_interface</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/boundary/cf_interface/</path>
    <filename>dir_2314a4444cbbca75b6e69be39f907534.html</filename>
    <file>CartCellDoubleLinearCFInterpolation.cpp</file>
    <file>CartCellDoubleQuadraticCFInterpolation.cpp</file>
    <file>CartSideDoubleQuadraticCFInterpolation.cpp</file>
    <file>CoarseFineBoundaryRefinePatchStrategy.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/coarsen_ops</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/coarsen_ops/</path>
    <filename>dir_3605c9ab889221fa673672650d36138a.html</filename>
    <file>CartCellDoubleCubicCoarsen.cpp</file>
    <file>CartSideDoubleCubicCoarsen.cpp</file>
    <file>LMarkerCoarsen.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/fortran</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/fortran/</path>
    <filename>dir_8b42ac1beded6aedffd3b5c5b79eaee3.html</filename>
    <file>minmod.f</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/include/ibtk</name>
    <path>/Users/griffith/code/IBAMR/ibtk/include/ibtk/</path>
    <filename>dir_8bf174eff2e89475adb76c315ddc2a3f.html</filename>
    <file>app_namespaces.h</file>
    <file>AppInitializer.h</file>
    <file>BGaussSeidelPreconditioner.h</file>
    <file>BJacobiPreconditioner.h</file>
    <file>CartCellDoubleBoundsPreservingConservativeLinearRefine.h</file>
    <file>CartCellDoubleCubicCoarsen.h</file>
    <file>CartCellDoubleLinearCFInterpolation.h</file>
    <file>CartCellDoubleQuadraticCFInterpolation.h</file>
    <file>CartCellDoubleQuadraticRefine.h</file>
    <file>CartCellRobinPhysBdryOp.h</file>
    <file>CartExtrapPhysBdryOp.h</file>
    <file>CartGridFunction.h</file>
    <file>CartGridFunctionSet.h</file>
    <file>CartSideDoubleCubicCoarsen.h</file>
    <file>CartSideDoubleDivPreservingRefine.h</file>
    <file>CartSideDoubleQuadraticCFInterpolation.h</file>
    <file>CartSideDoubleSpecializedConstantRefine.h</file>
    <file>CartSideDoubleSpecializedLinearRefine.h</file>
    <file>CartSideRobinPhysBdryOp.h</file>
    <file>CCLaplaceOperator.h</file>
    <file>CCPoissonHypreLevelSolver.h</file>
    <file>CCPoissonPETScLevelSolver.h</file>
    <file>CCPoissonPointRelaxationFACOperator.h</file>
    <file>CCPoissonSolverManager.h</file>
    <file>CellNoCornersFillPattern.h</file>
    <file>CoarseFineBoundaryRefinePatchStrategy.h</file>
    <file>CoarsenPatchStrategySet.h</file>
    <file>compiler_hints.h</file>
    <file>CopyToRootSchedule.h</file>
    <file>CopyToRootTransaction.h</file>
    <file>DebuggingUtilities.h</file>
    <file>EdgeDataSynchronization.h</file>
    <file>EdgeSynchCopyFillPattern.h</file>
    <file>ExtendedRobinBcCoefStrategy.h</file>
    <file>FaceDataSynchronization.h</file>
    <file>FaceSynchCopyFillPattern.h</file>
    <file>FACPreconditioner.h</file>
    <file>FACPreconditionerStrategy.h</file>
    <file>FEDataManager.h</file>
    <file>FixedSizedStream-inl.h</file>
    <file>FixedSizedStream.h</file>
    <file>GeneralOperator.h</file>
    <file>GeneralSolver.h</file>
    <file>HierarchyGhostCellInterpolation.h</file>
    <file>HierarchyIntegrator.h</file>
    <file>HierarchyMathOps.h</file>
    <file>ibtk.h</file>
    <file>IBTK_CHKERRQ.h</file>
    <file>ibtk_enums.h</file>
    <file>ibtk_utilities.h</file>
    <file>IndexUtilities-inl.h</file>
    <file>IndexUtilities.h</file>
    <file>JacobianOperator.h</file>
    <file>KrylovLinearSolver.h</file>
    <file>KrylovLinearSolverManager.h</file>
    <file>KrylovLinearSolverPoissonSolverInterface.h</file>
    <file>LaplaceOperator.h</file>
    <file>LData-inl.h</file>
    <file>LData.h</file>
    <file>LDataManager-inl.h</file>
    <file>LDataManager.h</file>
    <file>LEInteractor.h</file>
    <file>libmesh_utilities.h</file>
    <file>LIndexSetData-inl.h</file>
    <file>LIndexSetData.h</file>
    <file>LIndexSetDataFactory.h</file>
    <file>LIndexSetVariable.h</file>
    <file>LinearOperator.h</file>
    <file>LinearSolver.h</file>
    <file>LInitStrategy.h</file>
    <file>LMarker-inl.h</file>
    <file>LMarker.h</file>
    <file>LMarkerCoarsen.h</file>
    <file>LMarkerRefine.h</file>
    <file>LMarkerSet.h</file>
    <file>LMarkerSetData.h</file>
    <file>LMarkerSetDataFactory.h</file>
    <file>LMarkerSetDataIterator.h</file>
    <file>LMarkerSetVariable.h</file>
    <file>LMarkerTransaction.h</file>
    <file>LMarkerUtilities.h</file>
    <file>LMesh-inl.h</file>
    <file>LMesh.h</file>
    <file>LNode-inl.h</file>
    <file>LNode.h</file>
    <file>LNodeIndex-inl.h</file>
    <file>LNodeIndex.h</file>
    <file>LNodeIndexSet.h</file>
    <file>LNodeIndexSetData.h</file>
    <file>LNodeIndexSetDataFactory.h</file>
    <file>LNodeIndexSetDataIterator.h</file>
    <file>LNodeIndexSetVariable.h</file>
    <file>LNodeIndexTransaction.h</file>
    <file>LNodeSet.h</file>
    <file>LNodeSetData.h</file>
    <file>LNodeSetDataFactory.h</file>
    <file>LNodeSetDataIterator.h</file>
    <file>LNodeSetVariable.h</file>
    <file>LNodeTransaction.h</file>
    <file>LSet-inl.h</file>
    <file>LSet.h</file>
    <file>LSetData-inl.h</file>
    <file>LSetData.h</file>
    <file>LSetDataFactory.h</file>
    <file>LSetDataIterator-inl.h</file>
    <file>LSetDataIterator.h</file>
    <file>LSetVariable.h</file>
    <file>LSiloDataWriter.h</file>
    <file>LTransaction.h</file>
    <file>muParserCartGridFunction.h</file>
    <file>muParserRobinBcCoefs.h</file>
    <file>namespaces.h</file>
    <file>NewtonKrylovSolver.h</file>
    <file>NewtonKrylovSolverManager.h</file>
    <file>NodeDataSynchronization.h</file>
    <file>NodeSynchCopyFillPattern.h</file>
    <file>NormOps.h</file>
    <file>ParallelEdgeMap.h</file>
    <file>ParallelMap.h</file>
    <file>ParallelSet.h</file>
    <file>PatchMathOps.h</file>
    <file>PETScKrylovLinearSolver.h</file>
    <file>PETScKrylovPoissonSolver.h</file>
    <file>PETScLevelSolver.h</file>
    <file>PETScMatLOWrapper.h</file>
    <file>PETScMatUtilities.h</file>
    <file>PETScMFFDJacobianOperator.h</file>
    <file>PETScMultiVec.h</file>
    <file>PETScNewtonKrylovSolver.h</file>
    <file>PETScPCLSWrapper.h</file>
    <file>PETScSAMRAIVectorReal-inl.h</file>
    <file>PETScSAMRAIVectorReal.h</file>
    <file>PETScSNESFunctionGOWrapper.h</file>
    <file>PETScSNESJacobianJOWrapper.h</file>
    <file>PETScVecUtilities.h</file>
    <file>PhysicalBoundaryUtilities.h</file>
    <file>PoissonFACPreconditioner.h</file>
    <file>PoissonFACPreconditionerStrategy.h</file>
    <file>PoissonSolver.h</file>
    <file>PoissonUtilities.h</file>
    <file>RefinePatchStrategySet.h</file>
    <file>RobinPhysBdryPatchStrategy.h</file>
    <file>SCLaplaceOperator.h</file>
    <file>SCPoissonHypreLevelSolver.h</file>
    <file>SCPoissonPETScLevelSolver.h</file>
    <file>SCPoissonPointRelaxationFACOperator.h</file>
    <file>SCPoissonSolverManager.h</file>
    <file>SideDataSynchronization.h</file>
    <file>SideNoCornersFillPattern.h</file>
    <file>SideSynchCopyFillPattern.h</file>
    <file>StaggeredPhysicalBoundaryHelper.h</file>
    <file>StandardTagAndInitStrategySet.h</file>
    <file>Streamable.h</file>
    <file>StreamableFactory.h</file>
    <file>StreamableManager-inl.h</file>
    <file>StreamableManager.h</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk</name>
    <path>/Users/griffith/code/IBAMR/ibtk/</path>
    <filename>dir_9541794473d30fc5726c7b7c76923b64.html</filename>
    <dir>/Users/griffith/code/IBAMR/ibtk/include</dir>
    <dir>/Users/griffith/code/IBAMR/ibtk/src</dir>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/solvers/impls</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/solvers/impls/</path>
    <filename>dir_4cc51fd2f15c3b15d8388f9c74ab7945.html</filename>
    <file>BGaussSeidelPreconditioner.cpp</file>
    <file>BJacobiPreconditioner.cpp</file>
    <file>CCLaplaceOperator.cpp</file>
    <file>CCPoissonHypreLevelSolver.cpp</file>
    <file>CCPoissonPETScLevelSolver.cpp</file>
    <file>CCPoissonPointRelaxationFACOperator.cpp</file>
    <file>CCPoissonSolverManager.cpp</file>
    <file>FACPreconditioner.cpp</file>
    <file>KrylovLinearSolverManager.cpp</file>
    <file>KrylovLinearSolverPoissonSolverInterface.cpp</file>
    <file>LaplaceOperator.cpp</file>
    <file>NewtonKrylovSolverManager.cpp</file>
    <file>PETScKrylovLinearSolver.cpp</file>
    <file>PETScKrylovPoissonSolver.cpp</file>
    <file>PETScLevelSolver.cpp</file>
    <file>PETScMFFDJacobianOperator.cpp</file>
    <file>PETScMultiVec.cpp</file>
    <file>PETScNewtonKrylovSolver.cpp</file>
    <file>PoissonFACPreconditioner.cpp</file>
    <file>PoissonFACPreconditionerStrategy.cpp</file>
    <file>PoissonSolver.cpp</file>
    <file>SCLaplaceOperator.cpp</file>
    <file>SCPoissonHypreLevelSolver.cpp</file>
    <file>SCPoissonPETScLevelSolver.cpp</file>
    <file>SCPoissonPointRelaxationFACOperator.cpp</file>
    <file>SCPoissonSolverManager.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/include</name>
    <path>/Users/griffith/code/IBAMR/ibtk/include/</path>
    <filename>dir_4c0b489fb1e3d12ac48de2c424aa440d.html</filename>
    <dir>/Users/griffith/code/IBAMR/ibtk/include/ibtk</dir>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/solvers/interfaces</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/solvers/interfaces/</path>
    <filename>dir_e672d4976ff9aa0562cee66a7320e0d4.html</filename>
    <file>FACPreconditionerStrategy.cpp</file>
    <file>GeneralOperator.cpp</file>
    <file>GeneralSolver.cpp</file>
    <file>JacobianOperator.cpp</file>
    <file>KrylovLinearSolver.cpp</file>
    <file>LinearOperator.cpp</file>
    <file>LinearSolver.cpp</file>
    <file>NewtonKrylovSolver.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/lagrangian</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/lagrangian/</path>
    <filename>dir_a84b77d7055aca310a4dbfa8f930cf87.html</filename>
    <file>FEDataManager.cpp</file>
    <file>LData.cpp</file>
    <file>LDataManager.cpp</file>
    <file>LEInteractor.cpp</file>
    <file>LIndexSetData.cpp</file>
    <file>LIndexSetDataFactory.cpp</file>
    <file>LIndexSetVariable.cpp</file>
    <file>LInitStrategy.cpp</file>
    <file>LMarker.cpp</file>
    <file>LMesh.cpp</file>
    <file>LNode.cpp</file>
    <file>LNodeIndex.cpp</file>
    <file>LSet.cpp</file>
    <file>LSetData.cpp</file>
    <file>LSetDataFactory.cpp</file>
    <file>LSetDataIterator.cpp</file>
    <file>LSetVariable.cpp</file>
    <file>LSiloDataWriter.cpp</file>
    <file>LTransaction.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/math</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/math/</path>
    <filename>dir_349da0a2e5827dd0eae35f4e22289c1a.html</filename>
    <file>HierarchyMathOps.cpp</file>
    <file>PatchMathOps.cpp</file>
    <file>PETScMatUtilities.cpp</file>
    <file>PETScVecUtilities.cpp</file>
    <file>PoissonUtilities.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/boundary/physical_boundary</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/boundary/physical_boundary/</path>
    <filename>dir_d724c3704602df60336c90ef92eaf7c4.html</filename>
    <file>CartCellRobinPhysBdryOp.cpp</file>
    <file>CartExtrapPhysBdryOp.cpp</file>
    <file>CartSideRobinPhysBdryOp.cpp</file>
    <file>ExtendedRobinBcCoefStrategy.cpp</file>
    <file>muParserRobinBcCoefs.cpp</file>
    <file>PhysicalBoundaryUtilities.cpp</file>
    <file>RobinPhysBdryPatchStrategy.cpp</file>
    <file>StaggeredPhysicalBoundaryHelper.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/refine_ops</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/refine_ops/</path>
    <filename>dir_c282f1d762e0c3056df65e8b6bfbb8da.html</filename>
    <file>CartCellDoubleBoundsPreservingConservativeLinearRefine.cpp</file>
    <file>CartCellDoubleQuadraticRefine.cpp</file>
    <file>CartSideDoubleDivPreservingRefine.cpp</file>
    <file>CartSideDoubleSpecializedConstantRefine.cpp</file>
    <file>CartSideDoubleSpecializedLinearRefine.cpp</file>
    <file>LMarkerRefine.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/solvers</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/solvers/</path>
    <filename>dir_d315453bc067053d0a86066a5a84164f.html</filename>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/solvers/impls</dir>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/solvers/interfaces</dir>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/solvers/wrappers</dir>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/</path>
    <filename>dir_d59845d153a1f4aa708dcab91a175469.html</filename>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/boundary</dir>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/coarsen_ops</dir>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/fortran</dir>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/lagrangian</dir>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/math</dir>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/refine_ops</dir>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/solvers</dir>
    <dir>/Users/griffith/code/IBAMR/ibtk/src/utilities</dir>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/utilities</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/utilities/</path>
    <filename>dir_6bc406b0640767d292dfb2cc9840ef20.html</filename>
    <file>AppInitializer.cpp</file>
    <file>CartGridFunction.cpp</file>
    <file>CartGridFunctionSet.cpp</file>
    <file>CellNoCornersFillPattern.cpp</file>
    <file>CoarsenPatchStrategySet.cpp</file>
    <file>CopyToRootSchedule.cpp</file>
    <file>CopyToRootTransaction.cpp</file>
    <file>DebuggingUtilities.cpp</file>
    <file>EdgeDataSynchronization.cpp</file>
    <file>EdgeSynchCopyFillPattern.cpp</file>
    <file>FaceDataSynchronization.cpp</file>
    <file>FaceSynchCopyFillPattern.cpp</file>
    <file>FixedSizedStream.cpp</file>
    <file>HierarchyIntegrator.cpp</file>
    <file>IndexUtilities.cpp</file>
    <file>LMarkerUtilities.cpp</file>
    <file>muParserCartGridFunction.cpp</file>
    <file>NodeDataSynchronization.cpp</file>
    <file>NodeSynchCopyFillPattern.cpp</file>
    <file>NormOps.cpp</file>
    <file>ParallelEdgeMap.cpp</file>
    <file>ParallelMap.cpp</file>
    <file>ParallelSet.cpp</file>
    <file>RefinePatchStrategySet.cpp</file>
    <file>SideDataSynchronization.cpp</file>
    <file>SideNoCornersFillPattern.cpp</file>
    <file>SideSynchCopyFillPattern.cpp</file>
    <file>StandardTagAndInitStrategySet.cpp</file>
    <file>Streamable.cpp</file>
    <file>StreamableFactory.cpp</file>
    <file>StreamableManager.cpp</file>
  </compound>
  <compound kind="dir">
    <name>/Users/griffith/code/IBAMR/ibtk/src/solvers/wrappers</name>
    <path>/Users/griffith/code/IBAMR/ibtk/src/solvers/wrappers/</path>
    <filename>dir_0849e0e45f63e45bb99052423f0e7d70.html</filename>
    <file>PETScMatLOWrapper.cpp</file>
    <file>PETScPCLSWrapper.cpp</file>
    <file>PETScSAMRAIVectorReal.cpp</file>
    <file>PETScSNESFunctionGOWrapper.cpp</file>
    <file>PETScSNESJacobianJOWrapper.cpp</file>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title></title>
    <filename>index</filename>
  </compound>
</tagfile>
