//
//  MBEContext.m
//  MetalKernel
//
//  Created by asd on 28/03/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//


#import "MBEContext.h"
#import <math.h>

@implementation MetalDevice
+(instancetype)newDevice: (id<MTLDevice>)device { // init a device w/ library
    MetalDevice*md=[[super alloc]init];
    md.device=device;
    md.library=[device newDefaultLibrary];
    return md;
}
@end

@implementation MBEContext

+ (instancetype)newContext {
    MBEContext *c=[[super alloc]init];
    
    c.devices=[[NSMutableArray alloc]init];
    for (id<MTLDevice>device in MTLCopyAllDevices()) // init all devices
        [c.devices addObject:[MetalDevice newDevice:device]];
    c.nDevices = [c.devices count];
    
    return c;
}

-(id<MTLBuffer>)getBuffer: (id<MTLDevice>) device data:(void*)data length:(NSInteger)len {
    return [device newBufferWithBytes:data  length:len options:MTLResourceStorageModeShared];
}

-(void)domainColor: (uint32*)colorVect w:(int)w h:(int)h {

    int size=w*h;
    NSMutableArray<id<MTLBuffer>>*colorBuffers=[NSMutableArray arrayWithCapacity:_nDevices];
    
    for (int d=0; d<_nDevices; d++) {
        MetalDevice*md=_devices[d];
        // prepare shader
        md.kernelFunction = [md.library newFunctionWithName:@"generateDomainColor"];
        md.pipeline = [md.device  newComputePipelineStateWithFunction:md.kernelFunction error:nil];
        md.commandBuffer = [[md.device newCommandQueue] commandBuffer];
        md.commandEncoder = [md.commandBuffer computeCommandEncoder];
        [md.commandEncoder setComputePipelineState:md.pipeline];
        
        colorBuffers[d]=[self getBuffer: md.device data:colorVect length:(size/_nDevices)*sizeof(*colorVect)];
        
        [md.commandEncoder setBuffer:colorBuffers[d] offset:0 atIndex:0]; // shader parameters:
        // colorBuffer, current device, # devices, w, h
        [md.commandEncoder setBytes:&d length:sizeof(d) atIndex:1];
        [md.commandEncoder setBytes:&_nDevices length:sizeof(_nDevices) atIndex:2];
        
        [md.commandEncoder setBytes:&w length:sizeof(w) atIndex:3];
        [md.commandEncoder setBytes:&h length:sizeof(h) atIndex:4];
        
        // sqrt thExeWidth -> real thread number
        auto sz = MTLSizeMake( sqrt(md.pipeline.threadExecutionWidth),  1,  1);
        [md.commandEncoder dispatchThreadgroups:sz threadsPerThreadgroup:sz];
        [md.commandEncoder endEncoding];
        
        [md.commandBuffer commit ];
        [md.commandBuffer waitUntilCompleted];
    }
    
    // merge buffers
    for (int d=0; d<_nDevices; d++) {
        uint szdev=size/_nDevices;
        uint fd=d*szdev, szBytes=(((d==_nDevices-1) ? size : (d+1)*szdev)-fd)*sizeof(uint32);
        
        memcpy(colorVect + fd, [colorBuffers[d] contents], szBytes);
    }
    
}

+(NSTimeInterval) timeIt: (void (^) (void))block {
    NSDate *start = [NSDate date];
    block();
    return -[start timeIntervalSinceNow] * 1000.;
}

@end
