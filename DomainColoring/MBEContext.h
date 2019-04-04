//
//  MBEContext.h
//  MetalKernel
//
//  Created by asd on 28/03/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#ifndef MBEContext_h
#define MBEContext_h

#import <Metal/Metal.h>

@interface MetalDevice : NSObject // single metal device object
@property id<MTLDevice> device;
@property id<MTLLibrary> library;
@property id<MTLCommandQueue> commandQueue;
@property id<MTLFunction> kernelFunction;
@property id<MTLCommandBuffer> commandBuffer;
@property id<MTLComputeCommandEncoder> commandEncoder;
@property id<MTLComputePipelineState> pipeline;

+(instancetype)newDevice: (id<MTLDevice>)device;
@end


@interface MBEContext : NSObject // contains all available devices
@property NSMutableArray<MetalDevice*> *devices;
@property NSUInteger nDevices;

+ (instancetype)newContext;
-(void)domainColor: (uint32_t*)colorVect w:(int)w h:(int)h;
+(NSTimeInterval) timeIt : (void (^) (void))block;

@end




#endif /* MBEContext_h */
