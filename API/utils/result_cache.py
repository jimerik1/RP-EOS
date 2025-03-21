"""
Cache system for storing and retrieving calculation results.
"""
import hashlib
import json
import pickle
import logging
from pathlib import Path
from typing import Dict, Any, Tuple, Optional, List, Union

logger = logging.getLogger(__name__)

class ResultsCache:
    """Cache for storing calculation results to avoid redundant calculations."""
    
    def __init__(self, cache_dir: Optional[str] = None, memory_size: int = 64):
        """
        Initialize the cache system.
        
        Args:
            cache_dir: Directory for disk cache (None for memory-only)
            memory_size: Maximum number of results to keep in memory
        """
        # Memory cache
        self.memory_cache: Dict[str, Any] = {}
        self.memory_size = memory_size
        
        # Disk cache
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir and not self.cache_dir.exists():
            self.cache_dir.mkdir(parents=True)
            
        logger.info(f"Initialized result cache: memory_size={memory_size}, "
                   f"disk_cache={'enabled' if cache_dir else 'disabled'}")
    
    def get_key(self, composition: List[Dict[str, Any]], variables: Dict[str, Any], 
               properties: List[str], options: Dict[str, Any]) -> str:
        """
        Generate a unique cache key based on input parameters.
        
        Args:
            composition: List of fluid components and fractions
            variables: Dictionary of calculation variables
            properties: List of properties to calculate
            options: Additional calculation options
            
        Returns:
            Cache key string
        """
        # Create a canonical representation of inputs
        cache_dict = {
            'composition': sorted([{'fluid': c['fluid'], 'fraction': c['fraction']} 
                                  for c in composition], key=lambda x: x['fluid']),
            'variables': {k: {'range': v['range'], 'resolution': v['resolution']} 
                         for k, v in variables.items()},
            'properties': sorted(properties),
            'options': {k: v for k, v in options.items() 
                      if k not in ['use_cache', 'max_workers']}
        }
        
        # Create hash from the serialized dict
        key_str = json.dumps(cache_dict, sort_keys=True)
        return hashlib.md5(key_str.encode()).hexdigest()
    
    def get(self, composition: List[Dict[str, Any]], variables: Dict[str, Any],
           properties: List[str], options: Dict[str, Any]) -> Optional[Tuple]:
        """
        Try to get results from cache.
        
        Args:
            composition: List of fluid components and fractions
            variables: Dictionary of calculation variables
            properties: List of properties to calculate
            options: Additional calculation options
            
        Returns:
            Cached result tuple or None if not found
        """
        key = self.get_key(composition, variables, properties, options)
        
        # First check memory cache
        if key in self.memory_cache:
            logger.info(f"Cache hit (memory) for key {key[:8]}...")
            return self.memory_cache[key]
        
        # Then check disk cache if available
        if self.cache_dir:
            cache_file = self.cache_dir / f"{key}.pkl"
            if cache_file.exists():
                try:
                    with open(cache_file, 'rb') as f:
                        result = pickle.load(f)
                        # Also update memory cache
                        self._add_to_memory_cache(key, result)
                        logger.info(f"Cache hit (disk) for key {key[:8]}...")
                        return result
                except Exception as e:
                    logger.warning(f"Error reading cache file: {e}")
        
        logger.debug(f"Cache miss for key {key[:8]}...")
        return None
    
    def set(self, composition: List[Dict[str, Any]], variables: Dict[str, Any],
           properties: List[str], options: Dict[str, Any], result: Tuple) -> None:
        """
        Add results to cache.
        
        Args:
            composition: List of fluid components and fractions
            variables: Dictionary of calculation variables
            properties: List of properties to calculate
            options: Additional calculation options
            result: Result tuple to cache
        """
        key = self.get_key(composition, variables, properties, options)
        
        # Update memory cache
        self._add_to_memory_cache(key, result)
        
        # Update disk cache if available
        if self.cache_dir:
            cache_file = self.cache_dir / f"{key}.pkl"
            try:
                with open(cache_file, 'wb') as f:
                    pickle.dump(result, f)
                logger.debug(f"Result saved to disk cache: {key[:8]}...")
            except Exception as e:
                logger.warning(f"Error writing cache file: {e}")
    
    def _add_to_memory_cache(self, key: str, result: Tuple) -> None:
        """
        Add to memory cache with LRU behavior.
        
        Args:
            key: Cache key
            result: Result to cache
        """
        # If cache is full, remove oldest entry
        if len(self.memory_cache) >= self.memory_size:
            oldest_key = next(iter(self.memory_cache))
            del self.memory_cache[oldest_key]
            logger.debug(f"Removed oldest cache entry: {oldest_key[:8]}...")
        
        self.memory_cache[key] = result
        logger.debug(f"Added result to memory cache: {key[:8]}...")