__constant sampler_t linearSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

kernel void backprojectPixelDriven2DCL(
	global image2d_t sino,
	global float* grid,
	int gridSizeX,
	int gridSizeY, 
	float pixelSizeX,
	float pixelSizeY,
	float maxSIdx,
	float deltaS,
	float maxThetaIdx,
	float deltaTheta
	) {
	
	// compute x, y from thread idx
	const unsigned int x = get_global_id(0);// x index
	const unsigned int y = get_global_id(1);// y index
	
	if (x >= gridSizeX || y >= gridSizeY) {
		return;
	}
	
    int idx = x + y*gridSizeX;
	grid[idx] = 0;

		
    //loop over all angles
	for(int thetaIdx = 0; thetaIdx < maxThetaIdx; ++thetaIdx)
	{
			float theta = thetaIdx * deltaTheta;
			float sinTheta = sin(theta);
			float cosTheta = cos(theta);
	
			float2 pos = {-gridSizeX * pixelSizeX / 2 + x * pixelSizeX, -gridSizeY * pixelSizeY / 2 + y * pixelSizeY};
			float s = pos.x * cosTheta + pos.y *sinTheta;
			float sIdx = (s + maxSIdx * deltaS / 2) / deltaS;
			if(sIdx > maxSIdx - 1 || sIdx < 0)
			    continue;
			float2 sinoIdx = {sIdx, thetaIdx};
		    float val = read_imagef(sino, linearSampler, sinoIdx).x;
			val = val * M_PI_F / maxThetaIdx;
			grid[idx] += val;
	
	}
	

	
	return;

	
}